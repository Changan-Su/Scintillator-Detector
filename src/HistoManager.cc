//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file analysis/AnaEx01/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4ParticleGun.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <sstream>

namespace {
std::mutex gOutputPathMutex;
bool gSharedOutputPathInitialized = false;
std::string gSharedOutputFileBasePath;
std::string gSharedOutputDir;
bool gMetadataSourceWritten = false;

std::string FormatDoubleForFolder(G4double value)
{
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(3) << value;
  std::string out = oss.str();

  while (!out.empty() && out.back() == '0') {
    out.pop_back();
  }
  if (!out.empty() && out.back() == '.') {
    out.pop_back();
  }
  if (out.empty()) {
    out = "0";
  }

  for (char& c : out) {
    if (c == '.') {
      c = 'p';
    } else if (c == '-') {
      c = 'm';
    }
  }
  return out;
}

std::string BuildTimestamp()
{
  const auto now = std::chrono::system_clock::now();
  const std::time_t nowTime = std::chrono::system_clock::to_time_t(now);
  std::tm localTm{};
#ifdef _WIN32
  localtime_s(&localTm, &nowTime);
#else
  localtime_r(&nowTime, &localTm);
#endif

  const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                      now.time_since_epoch())
                      .count() %
                  1000;

  std::ostringstream ts;
  ts << std::put_time(&localTm, "%Y%m%d_%H%M%S") << "_" << std::setw(3)
     << std::setfill('0') << ms;
  return ts.str();
}

}  // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.



HistoManager::HistoManager()
{
  // Create or get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("csv");
  // the default file type can be overriden in run macro
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  if (!fFactoryOn) {
    analysisManager->SetVerboseLevel(1);

#ifdef G4MULTITHREADED
    analysisManager->SetNtupleMerging(true);
#endif

    analysisManager->SetHistoDirectoryName("histo");
    analysisManager->SetNtupleDirectoryName("ntuple");
  }

  if (!fOutputPathInitialized) {
    std::lock_guard<std::mutex> lock(gOutputPathMutex);
    if (!gSharedOutputPathInitialized) {
      const auto* detConstruction = static_cast<const B1::DetectorConstruction*>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());

      // Build folder name: [prefix_]<timestamp>
      std::string ts = BuildTimestamp();
      std::string prefix = detConstruction->GetResultsPrefix();
      std::string folderName = prefix.empty() ? ts : (prefix + "_" + ts);

      std::filesystem::path outputDir = std::filesystem::path("Results") / folderName;
      if (std::filesystem::exists(outputDir)) {
        int dupIdx = 1;
        while (std::filesystem::exists(outputDir)) {
          outputDir = std::filesystem::path("Results")
                    / (folderName + "_dup" + std::to_string(dupIdx));
          ++dupIdx;
        }
      }
      std::filesystem::create_directories(outputDir);

      // Write metadata.csv with detector parameters (source params appended later by worker thread)
      std::filesystem::path metaPath = outputDir / "metadata.csv";
      std::ofstream meta(metaPath.string());
      if (meta.is_open()) {
        meta << "Key,Value\n";
        // detector parameters
        meta << "detector.surfaceSigma," << detConstruction->GetSurfaceSigma() << "\n";
        meta << "detector.arrayNx," << detConstruction->GetCrystal_nx() << "\n";
        meta << "detector.arrayNy," << detConstruction->GetCrystal_ny() << "\n";
        meta << "detector.arrayNz," << detConstruction->GetCrystal_nz() << "\n";
        meta << "detector.crystalGap_mm," << detConstruction->GetCrystal_gap() / mm << "\n";
        meta << "detector.crystalSize_mm," << detConstruction->Getcrystal_l() / mm << "\n";
        meta << "detector.crystalSizeY_mm," << detConstruction->Getcrystal_ly() / mm << "\n";
        meta << "detector.fillterRatioY," << detConstruction->GetFillterRatioY() << "\n";
        meta << "detector.fillterRatioZ," << detConstruction->GetFillterRatioZ() << "\n";
        meta << "detector.fillterPosRatioY," << detConstruction->GetFillterPosRatioY() << "\n";
        meta << "detector.fillterPosRatioZ," << detConstruction->GetFillterPosRatioZ() << "\n";
        meta << "results.prefix," << prefix << "\n";
        // source.* params are appended by the first worker thread (genAction is null on master)
        meta.close();
        G4cout << "\n----> metadata.csv written to " << metaPath.string() << G4endl;
      }

      gSharedOutputFileBasePath = (outputDir / "AnaEx01").string();
      gSharedOutputDir = outputDir.string();
      gSharedOutputPathInitialized = true;
    }

    fOutputFileBasePath = gSharedOutputFileBasePath;
    fOutputPathInitialized = true;
  }

  // Append source params to metadata once a worker thread is available
  // (PrimaryGeneratorAction is null on the master thread)
  {
    std::lock_guard<std::mutex> lock(gOutputPathMutex);
    if (!gMetadataSourceWritten) {
      const auto* genAction = static_cast<const B1::PrimaryGeneratorAction*>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
      if (genAction && !gSharedOutputDir.empty()) {
        std::filesystem::path metaPath = std::filesystem::path(gSharedOutputDir) / "metadata.csv";
        std::ofstream meta(metaPath.string(), std::ios::app);
        if (meta.is_open()) {
          meta << "source.mode," << genAction->GetSourceMode() << "\n";
          meta << "source.distribution," << genAction->GetSourceDistribution() << "\n";
          G4ThreeVector fp = genAction->GetFpSource();
          meta << "source.fp_source_x_mm," << fp.x() / mm << "\n";
          meta << "source.fp_source_y_mm," << fp.y() / mm << "\n";
          meta << "source.fp_source_z_mm," << fp.z() / mm << "\n";
          meta.close();
        }
        gMetadataSourceWritten = true;
      }
    }
  }

  G4bool fileOpen = analysisManager->OpenFile(fOutputFileBasePath);
  if (!fileOpen) {
    G4cerr << "\n---> HistoManager::Book(): cannot open " << analysisManager->GetFileName()
           << G4endl;
    return;
  }

  if (!fFactoryOn) {
    // Create histograms
    analysisManager->CreateH1("EAbs", "Edep in absorber (MeV)", 100, 0., 800 * MeV);
    analysisManager->CreateH1("EGap", "Edep in gap (MeV)", 100, 0., 100 * MeV);
    analysisManager->CreateH1("LAbs", "trackL in absorber (mm)", 100, 0., 1 * m);
    analysisManager->CreateH1("LGap", "trackL in gap (mm)", 100, 0., 50 * cm);
    analysisManager->CreateH1("PhotonCounts", "Number of optical photons in PMT", 200, 0, 2000);
    analysisManager->CreateH1("PhotonDepth","Depth index of captured photons", 20, -0.5, 19.5);
    // Create Ntuples
    analysisManager->CreateNtuple("Ntuple1", "Edep");
    analysisManager->CreateNtupleDColumn("Eabs");  // column Id = 0
    analysisManager->CreateNtupleDColumn("Egap");  // column Id = 1
    analysisManager->FinishNtuple();

    analysisManager->CreateNtuple("Ntuple2", "TrackL");
    analysisManager->CreateNtupleDColumn("Labs");  // column Id = 0
    analysisManager->CreateNtupleDColumn("Lgap");  // column Id = 1
    analysisManager->FinishNtuple();

    // Legacy photon-related CSV outputs are temporarily disabled.
    // Kept here as commented code for easy rollback if needed.
    //
    // analysisManager->CreateNtuple("PhotonNtuple", "PhotonCounts");
    // analysisManager->CreateNtupleIColumn("PhotonCounts");
    // analysisManager->FinishNtuple();
    //
    // analysisManager->CreateNtuple("PhotonGenerated", "GeneratedCounts");
    // analysisManager->CreateNtupleIColumn("GeneratedCounts");
    // analysisManager->FinishNtuple();
    //
    // fDepthNtupleId = analysisManager->CreateNtuple("PhotonDepthNtuple", "Depth index of each captured photon");
    // analysisManager->CreateNtupleIColumn("DepthIndex");
    // analysisManager->FinishNtuple();
    //
    // analysisManager->CreateNtuple("PhotonLeft", "PhotonCounts_Left");
    // analysisManager->CreateNtupleIColumn("PhotonCounts_Left");
    // analysisManager->FinishNtuple();
    //
    // analysisManager->CreateNtuple("PhotonRight", "PhotonCounts_Right");
    // analysisManager->CreateNtupleIColumn("PhotonCounts_Right");
    // analysisManager->FinishNtuple();
    //
    // fPhotonLRNtupleId = analysisManager->CreateNtuple("PhotonLRPerRod",
    //                                               "Per-rod Left/Right photon counts per event");
    // analysisManager->CreateNtupleIColumn("EventID");
    // analysisManager->CreateNtupleIColumn("iz");
    // analysisManager->CreateNtupleIColumn("iy");
    // analysisManager->CreateNtupleIColumn("Left");
    // analysisManager->CreateNtupleIColumn("Right");
    // analysisManager->FinishNtuple();

    fPhotonFaceEventNtupleId = analysisManager->CreateNtuple(
        "PhotonFaceBlockEvent", "Per-event photon counts on each SiPM block");
    analysisManager->CreateNtupleIColumn("EventID");
    analysisManager->CreateNtupleIColumn("CrystalID");
    analysisManager->CreateNtupleIColumn("iy");
    analysisManager->CreateNtupleIColumn("iz");
    analysisManager->CreateNtupleIColumn("Face");
    analysisManager->CreateNtupleIColumn("j");
    analysisManager->CreateNtupleIColumn("k");
    analysisManager->CreateNtupleIColumn("SiPMBlockID");
    analysisManager->CreateNtupleIColumn("PhotonCount");
    analysisManager->FinishNtuple();

    fFactoryOn = true;
  }

  G4cout << "\n----> Output file is open in " << analysisManager->GetFileName() << "."
         << analysisManager->GetFileType() << G4endl;



  // G4cout << "PhotonDepthNtuple rows: " 
  //      << G4AnalysisManager::Instance()->GetNtuple(4)->GetNtupleRowCount()
  //      << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




void HistoManager::FillPhotonNtuple(G4int photonCounts)
{
    // Legacy output path is disabled.
    (void)photonCounts;
    return;
}

void HistoManager::FillPhotonGeneratedNtuple(G4int photonGenerated)
{
    // Legacy output path is disabled.
    (void)photonGenerated;
    return;
}

void HistoManager::FillPhotonHisto(G4int photonCounts)
  {G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillH1(4, photonCounts); }

void HistoManager::Save()
{
  if (!fFactoryOn) {
    return;
  }

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  analysisManager->Write();
  analysisManager->CloseFile();

  G4cout << "\n----> Histograms and ntuples are saved\n" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillHisto(G4int ih, G4double xbin, G4double weight)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillH1(ih, xbin, weight);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Normalize(G4int ih, G4double fac)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  auto h1 = analysisManager->GetH1(ih);
  if (h1 != nullptr) {
    h1->scale(fac);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void HistoManager::FillNtuple(G4double energyAbs, G4double energyGap, G4double trackLAbs,G4double trackLGap)
// {
//   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
//   // Fill 1st ntuple ( id = 0)
//   analysisManager->FillNtupleDColumn(0, 0, energyAbs);
//   analysisManager->FillNtupleDColumn(0, 1, energyGap);
//   analysisManager->AddNtupleRow(0);
//   // Fill 2nd ntuple ( id = 1)
//   analysisManager->FillNtupleDColumn(1, 0, trackLAbs);
//   analysisManager->FillNtupleDColumn(1, 1, trackLGap);
//   analysisManager->AddNtupleRow(1);
// }



void HistoManager::FillNtuple(G4double energyAbs, G4double trackLAbs)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Fill energy (Ntuple1)
  analysisManager->FillNtupleDColumn(0, 0, energyAbs);
  analysisManager->FillNtupleDColumn(0, 1, 0.);  // 默认 gap 为 0
  analysisManager->AddNtupleRow(0);

  // Fill track length (Ntuple2)
  analysisManager->FillNtupleDColumn(1, 0, trackLAbs);
  analysisManager->FillNtupleDColumn(1, 1, 0.);  // 默认 gap 为 0
  analysisManager->AddNtupleRow(1);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void HistoManager::PrintStatistic()
{

  if (!fFactoryOn) {
    return;
  }

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4cout << "\n ----> print histograms statistic \n" << G4endl;
  for (G4int i = 0; i < analysisManager->GetNofH1s(); ++i) {
    G4String name = analysisManager->GetH1Name(i);
    auto h1 = analysisManager->GetH1(i);

    G4String unitCategory;
    if (name[0U] == 'E') {
      unitCategory = "Energy";
    }
    if (name[0U] == 'L') {
      unitCategory = "Length";
    }
    // we use an explicit unsigned int type for operator [] argument
    // to avoid problems with windows compiler

    G4cout << name << ": mean = " << G4BestUnit(h1->mean(), unitCategory)
           << " rms = " << G4BestUnit(h1->rms(), unitCategory) << G4endl;
  }
}

void HistoManager::FillPhotonDepth(G4int depthIdx)
{   
    if (fDepthNtupleId < 0) {
      return;
    }
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(fDepthNtupleId, 0, depthIdx);  // 填写 depth
    analysisManager->AddNtupleRow(fDepthNtupleId);                    // 添加到 ntuple       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillPhotonLeft(G4int photonCounts)
{
    (void)photonCounts;
    // Legacy output path is disabled.
}

void HistoManager::FillPhotonRight(G4int photonCounts)
{
    (void)photonCounts;
    // Legacy output path is disabled.
}


void HistoManager::FillPhotonLRPerRod(G4int iz, G4int iy,
                                      G4int left, G4int right,
                                      G4int eventId)
{
    (void)iz;
    (void)iy;
    (void)left;
    (void)right;
    (void)eventId;
    // Legacy output path is disabled.
}

void HistoManager::FillPhotonFaceBlockEvent(G4int eventId, G4int crystalId, G4int iy, G4int iz,
                                            G4int face, G4int j, G4int k, G4int sipmBlockId,
                                            G4int photonCount)
{
    if (photonCount <= 0) {
      return;
    }

    auto* am = G4AnalysisManager::Instance();
    if (fPhotonFaceEventNtupleId >= 0) {
      am->FillNtupleIColumn(fPhotonFaceEventNtupleId, 0, eventId);
      am->FillNtupleIColumn(fPhotonFaceEventNtupleId, 1, crystalId);
      am->FillNtupleIColumn(fPhotonFaceEventNtupleId, 2, iy);
      am->FillNtupleIColumn(fPhotonFaceEventNtupleId, 3, iz);
      am->FillNtupleIColumn(fPhotonFaceEventNtupleId, 4, face);
      am->FillNtupleIColumn(fPhotonFaceEventNtupleId, 5, j);
      am->FillNtupleIColumn(fPhotonFaceEventNtupleId, 6, k);
      am->FillNtupleIColumn(fPhotonFaceEventNtupleId, 7, sipmBlockId);
      am->FillNtupleIColumn(fPhotonFaceEventNtupleId, 8, photonCount);
      am->AddNtupleRow(fPhotonFaceEventNtupleId);
    }

    // Block total CSV has been removed by request.
}

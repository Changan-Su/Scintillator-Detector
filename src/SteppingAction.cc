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
//
/// \file B1/src/SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"
#include "G4OpticalPhoton.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"


#include "G4Event.hh"
#pragma message ("!!! Using EventAction.hh from: " __FILE__)
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction) : fEventAction(eventAction) {}


void SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume || !flogicSiPM) {
    const auto detConstruction = static_cast<const DetectorConstruction*>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detConstruction->GetScoringVolume();
    flogicSiPM = detConstruction->GetlogicSiPM();
  }

  G4Track* track = step->GetTrack();


  if (track->GetDefinition() == G4OpticalPhoton::Definition()) {
    // 只在刚生成的那一刻+1
    if (track->GetCurrentStepNumber() == 1) {
        fEventAction->AddPhotonGenerated(1);
    }
    G4LogicalVolume* volume =
        step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    G4String volName = volume->GetName();
    // G4cout << "[DEBUG] Optical photon at volume: " << volName << G4endl;
    // if (volName.contains("SiPM")) {
    //     G4int copyNo_SiPM = step->GetPreStepPoint()
    //                             ->GetTouchableHandle()->GetCopyNumber();

    //     G4int DOI_iz = copyNo_SiPM / 1000000;
    //     fEventAction->AddPhotonDepth(DOI_iz);

    //     fEventAction->AddPhotonInPMT(1);//Photons Counts

    //     if (copyNo_SiPM % 100 == 1) {
    //     // Left SiPM
    //     // G4int FuckU = 1;
    //     fEventAction->AddPhotonLeft(1);
    //     } 
    //     else if (copyNo_SiPM % 100 == 2) {
    //     // Right SiPM
    //     fEventAction->AddPhotonRight(1);
    //     }

    //     // G4cout<< "Debug:DOI" << DOI_iz << G4endl;
    //     track->SetTrackStatus(fStopAndKill);  // 防止重复计数
    // }

    if (volName.contains("SiPM")) {
    const G4int copyNo = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
    const G4int iz = copyNo / 1000000;
    const G4int iy = (copyNo - iz*1000000) / 100;   // 提取 iy
    const G4int endTag = copyNo % 100;          // 1=Left, 2=Right

    if (endTag == 1) {
        fEventAction->AddPhotonLeftAt(iy, iz);
    } else if (endTag == 2) {
        fEventAction->AddPhotonRightAt(iy, iz);
    }
    track->SetTrackStatus(fStopAndKill);
}

    // G4ThreeVector pos = track->GetPosition();
    // // 你可以根据晶体阵列的实际大小设置阈值，比如 ±Crystal_x/2 ±1cm 边界
    // G4double maxX = fDetectorConstruction->GetCrystal_x()/2 + 0.1 * fDetectorConstruction->Getcrystal_l();
    // G4double maxY = fDetectorConstruction->GetCrystal_y()/2;
    // G4double maxZ = fDetectorConstruction->GetCrystal_z()/2;

    // if (std::abs(pos.x()) > maxX ||
    //     std::abs(pos.y()) > maxY ||
    //     std::abs(pos.z()) > maxZ) {
    //     track->SetTrackStatus(fStopAndKill);  // 超出范围直接杀
    // }
}

  G4double edepStep   = step->GetTotalEnergyDeposit();
  G4double stepLength = step->GetStepLength();

  if (edepStep > 0.) {
    fEventAction->AddEdep(edepStep);
    fEventAction->AddAbsorption(edepStep, stepLength);
    //  G4cout << "Edep: " << edepStep / keV << " keV in volume: " 
    //        << step->GetPreStepPoint()->GetPhysicalVolume()->GetName() << G4endl;
      const auto* secondaries = step->GetSecondaryInCurrentStep();
    int nOptical = 0;
    for (const auto& s : *secondaries) {
        if (s->GetDefinition() == G4OpticalPhoton::Definition()) {
            nOptical++;
        }
    }
    // G4cout << "[ScintDebug] Edep: " << edepStep / keV 
    //        << " keV  => OpticalPhotons: " << nOptical << G4endl;
  
  
          }
}  

}  // namespace B1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


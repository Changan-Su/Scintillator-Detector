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
/// \file B1/src/EventAction.cc
/// \brief Implementation of the B1::EventAction class

#include "EventAction.hh"
#include "G4Event.hh"
#include "RunAction.hh"
#include "HistoManager.hh"
#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"

#include "Randomize.hh"//随机数库
#include "DetectorConstruction.hh"

namespace B1
{

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// EventAction::EventAction(RunAction* runAction) : fRunAction(runAction) {}

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void EventAction::BeginOfEventAction(const G4Event*)
// {
//   fEdep = 0.;
// }


// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void EventAction::EndOfEventAction(const G4Event*)
// {
//   // accumulate statistics in run action
//   fRunAction->AddEdep(fEdep);
// }

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



EventAction::EventAction(RunAction* runAction,HistoManager* histo) : fRunAction(runAction),fHistoManager(histo){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  fEdep = 0.;
  fEnergyAbs  = 0.;
  fTrackLAbs  = 0.;
  fPhotonCountDetector = 0;
  fPhotonCountPMT = 0;
  fPhotonCountGenerated = 0;
  fPhotonCountLeft = 0;
  fPhotonCountRight = 0;


  G4int evtNb = evt->GetEventID();
  if (evtNb % fPrintModulo == 0) {
    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
  }

    const auto det = static_cast<const DetectorConstruction*>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fNy = det->GetCrystal_ny();
    fNz = det->GetCrystal_nz();
    fLeftPerRod.assign(fNy*fNz, 0);
    fRightPerRod.assign(fNy*fNz, 0);

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{

    G4double edepToRecord = fEdep;
  
    if (fApplySmearing && fEdep > 0.) {
        G4double fwhm = 0.07;  // 7%分辨率
        
        G4double sigma = fEdep * fwhm / 2.355;
        edepToRecord = CLHEP::RandGaussQ::shoot(fEdep, sigma);
    }


    for (auto d : fPhotonDOIs) 
    {
    fHistoManager->FillPhotonDepth(d);
    }
    fPhotonDOIs.clear();


    fRunAction->AddEdep(fEdep);
    fRunAction->FillPerEvent(fEnergyAbs /MeV, fTrackLAbs /mm);

    fHistoManager->FillHisto(0, edepToRecord);          // 展宽or不展宽
    fHistoManager->FillHisto(1, fTrackLAbs);
    fHistoManager->FillNtuple(edepToRecord, fTrackLAbs); 
    //Photon统计
    G4int photonToRecord = fPhotonCountPMT; 

    fHistoManager->FillPhotonHisto(photonToRecord); 
    fHistoManager->FillPhotonNtuple(fPhotonCountPMT);
    fHistoManager->FillPhotonGeneratedNtuple(fPhotonCountGenerated);
    fHistoManager->FillPhotonLeft(fPhotonCountLeft);
    fHistoManager->FillPhotonRight(fPhotonCountRight);

    const int eventId = evt->GetEventID();
    for (int iz=0; iz<fNz; ++iz)
      for (int iy=0; iy<fNy; ++iy) {
        int L = fLeftPerRod[RodIndex(iy,iz)];
        int R = fRightPerRod[RodIndex(iy,iz)];
        fHistoManager->FillPhotonLRPerRod(iz, iy, L, R, eventId);
      }

  }


}  // namespace B1

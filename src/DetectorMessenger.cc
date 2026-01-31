//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// ********************************************************************
//
/// \file DetectorMessenger.cc
/// \brief Implementation of the B1::DetectorMessenger class

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

namespace B1
{

DetectorMessenger::DetectorMessenger(DetectorConstruction* det)
  : G4UImessenger(), fDetector(det)
{
  fDetDir = new G4UIdirectory("/detector/");
  fDetDir->SetGuidance("Crystal array and gap (geometry) commands.");

  fArrayNxCmd = new G4UIcmdWithAnInteger("/detector/arrayNx", this);
  fArrayNxCmd->SetGuidance("Number of crystals in x direction.");
  fArrayNxCmd->SetParameterName("Nx", false);
  fArrayNxCmd->SetRange("Nx >= 1");
  fArrayNxCmd->SetDefaultValue(11);

  fArrayNyCmd = new G4UIcmdWithAnInteger("/detector/arrayNy", this);
  fArrayNyCmd->SetGuidance("Number of crystals in y direction.");
  fArrayNyCmd->SetParameterName("Ny", false);
  fArrayNyCmd->SetRange("Ny >= 1");
  fArrayNyCmd->SetDefaultValue(7);

  fArrayNzCmd = new G4UIcmdWithAnInteger("/detector/arrayNz", this);
  fArrayNzCmd->SetGuidance("Number of crystals in z direction.");
  fArrayNzCmd->SetParameterName("Nz", false);
  fArrayNzCmd->SetRange("Nz >= 1");
  fArrayNzCmd->SetDefaultValue(7);

  fCrystalGapCmd = new G4UIcmdWithADoubleAndUnit("/detector/crystalGap", this);
  fCrystalGapCmd->SetGuidance("Gap between crystals.");
  fCrystalGapCmd->SetParameterName("Gap", false);
  fCrystalGapCmd->SetRange("Gap >= 0.");
  fCrystalGapCmd->SetDefaultValue(0.1);
  fCrystalGapCmd->SetDefaultUnit("mm");

  fCrystalSizeCmd = new G4UIcmdWithADoubleAndUnit("/detector/crystalSize", this);
  fCrystalSizeCmd->SetGuidance("Single crystal size (x/z dimension).");
  fCrystalSizeCmd->SetParameterName("Size", false);
  fCrystalSizeCmd->SetRange("Size > 0.");
  fCrystalSizeCmd->SetDefaultValue(3.);
  fCrystalSizeCmd->SetDefaultUnit("mm");

  fCrystalSizeYCmd = new G4UIcmdWithADoubleAndUnit("/detector/crystalSizeY", this);
  fCrystalSizeYCmd->SetGuidance("Single crystal size in y dimension.");
  fCrystalSizeYCmd->SetParameterName("SizeY", false);
  fCrystalSizeYCmd->SetRange("SizeY > 0.");
  fCrystalSizeYCmd->SetDefaultValue(3.);
  fCrystalSizeYCmd->SetDefaultUnit("mm");

  fFillterRatioYCmd = new G4UIcmdWithADouble("/detector/fillterRatioY", this);
  fFillterRatioYCmd->SetGuidance("Fillter size ratio in Y direction (0-1, relative to crystal_ly).");
  fFillterRatioYCmd->SetParameterName("RatioY", false);
  fFillterRatioYCmd->SetRange("RatioY >= 0. && RatioY <= 1.");
  fFillterRatioYCmd->SetDefaultValue(0.3);

  fFillterRatioZCmd = new G4UIcmdWithADouble("/detector/fillterRatioZ", this);
  fFillterRatioZCmd->SetGuidance("Fillter size ratio in Z direction (0-1, relative to crystal_l).");
  fFillterRatioZCmd->SetParameterName("RatioZ", false);
  fFillterRatioZCmd->SetRange("RatioZ >= 0. && RatioZ <= 1.");
  fFillterRatioZCmd->SetDefaultValue(1.0);

  fFillterPosRatioYCmd = new G4UIcmdWithADouble("/detector/fillterPosRatioY", this);
  fFillterPosRatioYCmd->SetGuidance("Fillter position ratio in Y (0-1, from center of crystal to center of fillter).");
  fFillterPosRatioYCmd->SetParameterName("PosRatioY", false);
  fFillterPosRatioYCmd->SetRange("PosRatioY >= 0. && PosRatioY <= 1.");
  fFillterPosRatioYCmd->SetDefaultValue(0.7);

  fFillterPosRatioZCmd = new G4UIcmdWithADouble("/detector/fillterPosRatioZ", this);
  fFillterPosRatioZCmd->SetGuidance("Fillter position ratio in Z (0-1, from center of crystal to center of fillter).");
  fFillterPosRatioZCmd->SetParameterName("PosRatioZ", false);
  fFillterPosRatioZCmd->SetRange("PosRatioZ >= 0. && PosRatioZ <= 1.");
  fFillterPosRatioZCmd->SetDefaultValue(0.0);

  fUpdateCmd = new G4UIcmdWithoutParameter("/detector/update", this);
  fUpdateCmd->SetGuidance("Update geometry: rebuild detector with current parameters.");
  fUpdateCmd->SetGuidance("Use this after changing array size, gap, or crystal dimensions.");
}

DetectorMessenger::~DetectorMessenger()
{
  delete fUpdateCmd;
  delete fFillterPosRatioZCmd;
  delete fFillterPosRatioYCmd;
  delete fFillterRatioZCmd;
  delete fFillterRatioYCmd;
  delete fArrayNxCmd;
  delete fArrayNyCmd;
  delete fArrayNzCmd;
  delete fCrystalGapCmd;
  delete fCrystalSizeCmd;
  delete fCrystalSizeYCmd;
  delete fDetDir;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fArrayNxCmd)
    fDetector->SetArrayNx(fArrayNxCmd->GetNewIntValue(newValue));
  else if (command == fArrayNyCmd)
    fDetector->SetArrayNy(fArrayNyCmd->GetNewIntValue(newValue));
  else if (command == fArrayNzCmd)
    fDetector->SetArrayNz(fArrayNzCmd->GetNewIntValue(newValue));
  else if (command == fCrystalGapCmd)
    fDetector->SetCrystalGap(fCrystalGapCmd->GetNewDoubleValue(newValue));
  else if (command == fCrystalSizeCmd)
    fDetector->SetCrystalSize(fCrystalSizeCmd->GetNewDoubleValue(newValue));
  else if (command == fCrystalSizeYCmd)
    fDetector->SetCrystalSizeY(fCrystalSizeYCmd->GetNewDoubleValue(newValue));
  else if (command == fFillterRatioYCmd)
    fDetector->SetFillterRatioY(fFillterRatioYCmd->GetNewDoubleValue(newValue));
  else if (command == fFillterRatioZCmd)
    fDetector->SetFillterRatioZ(fFillterRatioZCmd->GetNewDoubleValue(newValue));
  else if (command == fFillterPosRatioYCmd)
    fDetector->SetFillterPosRatioY(fFillterPosRatioYCmd->GetNewDoubleValue(newValue));
  else if (command == fFillterPosRatioZCmd)
    fDetector->SetFillterPosRatioZ(fFillterPosRatioZCmd->GetNewDoubleValue(newValue));
  else if (command == fUpdateCmd) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
    G4cout << "Geometry updated with current parameters." << G4endl;
  }
}

}  // namespace B1

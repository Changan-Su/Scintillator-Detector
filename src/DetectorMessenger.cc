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
#include "G4UIcmdWithABool.hh"
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

  fUseGeometryMacCmd = new G4UIcmdWithABool("/detector/useGeometryMac", this);
  fUseGeometryMacCmd->SetGuidance("If true, crystal array uses geometry from geometry.mac (UI commands). If false, uses manual parameters.");
  fUseGeometryMacCmd->SetParameterName("useMac", false);
  fUseGeometryMacCmd->SetDefaultValue(true);

  fManualArrayNxCmd = new G4UIcmdWithAnInteger("/detector/manualArrayNx", this);
  fManualArrayNxCmd->SetGuidance("Manual: number of crystals in x (used when useGeometryMac false).");
  fManualArrayNxCmd->SetParameterName("Nx", false);
  fManualArrayNxCmd->SetRange("Nx >= 1");
  fManualArrayNxCmd->SetDefaultValue(9);
  fManualArrayNyCmd = new G4UIcmdWithAnInteger("/detector/manualArrayNy", this);
  fManualArrayNyCmd->SetGuidance("Manual: number of crystals in y.");
  fManualArrayNyCmd->SetParameterName("Ny", false);
  fManualArrayNyCmd->SetRange("Ny >= 1");
  fManualArrayNyCmd->SetDefaultValue(1);
  fManualArrayNzCmd = new G4UIcmdWithAnInteger("/detector/manualArrayNz", this);
  fManualArrayNzCmd->SetGuidance("Manual: number of crystals in z.");
  fManualArrayNzCmd->SetParameterName("Nz", false);
  fManualArrayNzCmd->SetRange("Nz >= 1");
  fManualArrayNzCmd->SetDefaultValue(1);
  fManualCrystalGapCmd = new G4UIcmdWithADoubleAndUnit("/detector/manualCrystalGap", this);
  fManualCrystalGapCmd->SetGuidance("Manual: gap between crystals (mm).");
  fManualCrystalGapCmd->SetParameterName("Gap", false);
  fManualCrystalGapCmd->SetRange("Gap >= 0.");
  fManualCrystalGapCmd->SetDefaultValue(0.);
  fManualCrystalGapCmd->SetDefaultUnit("mm");
  fManualCrystalSizeCmd = new G4UIcmdWithADoubleAndUnit("/detector/manualCrystalSize", this);
  fManualCrystalSizeCmd->SetGuidance("Manual: crystal size x/z (mm).");
  fManualCrystalSizeCmd->SetParameterName("Size", false);
  fManualCrystalSizeCmd->SetRange("Size > 0.");
  fManualCrystalSizeCmd->SetDefaultValue(3.);
  fManualCrystalSizeCmd->SetDefaultUnit("mm");
  fManualCrystalSizeYCmd = new G4UIcmdWithADoubleAndUnit("/detector/manualCrystalSizeY", this);
  fManualCrystalSizeYCmd->SetGuidance("Manual: crystal size y (mm).");
  fManualCrystalSizeYCmd->SetParameterName("SizeY", false);
  fManualCrystalSizeYCmd->SetRange("SizeY > 0.");
  fManualCrystalSizeYCmd->SetDefaultValue(3.);
  fManualCrystalSizeYCmd->SetDefaultUnit("mm");
  fManualFillterRatioYCmd = new G4UIcmdWithADouble("/detector/manualFillterRatioY", this);
  fManualFillterRatioYCmd->SetGuidance("Manual: fillter size ratio Y (0-1).");
  fManualFillterRatioYCmd->SetParameterName("RatioY", false);
  fManualFillterRatioYCmd->SetRange("RatioY >= 0. && RatioY <= 1.");
  fManualFillterRatioYCmd->SetDefaultValue(0.3);
  fManualFillterRatioZCmd = new G4UIcmdWithADouble("/detector/manualFillterRatioZ", this);
  fManualFillterRatioZCmd->SetGuidance("Manual: fillter size ratio Z (0-1).");
  fManualFillterRatioZCmd->SetParameterName("RatioZ", false);
  fManualFillterRatioZCmd->SetRange("RatioZ >= 0. && RatioZ <= 1.");
  fManualFillterRatioZCmd->SetDefaultValue(1.0);
  fManualFillterPosRatioYCmd = new G4UIcmdWithADouble("/detector/manualFillterPosRatioY", this);
  fManualFillterPosRatioYCmd->SetGuidance("Manual: fillter position ratio Y (0-1).");
  fManualFillterPosRatioYCmd->SetParameterName("PosRatioY", false);
  fManualFillterPosRatioYCmd->SetRange("PosRatioY >= 0. && PosRatioY <= 1.");
  fManualFillterPosRatioYCmd->SetDefaultValue(0.7);
  fManualFillterPosRatioZCmd = new G4UIcmdWithADouble("/detector/manualFillterPosRatioZ", this);
  fManualFillterPosRatioZCmd->SetGuidance("Manual: fillter position ratio Z (0-1).");
  fManualFillterPosRatioZCmd->SetParameterName("PosRatioZ", false);
  fManualFillterPosRatioZCmd->SetRange("PosRatioZ >= 0. && PosRatioZ <= 1.");
  fManualFillterPosRatioZCmd->SetDefaultValue(0.0);
}

DetectorMessenger::~DetectorMessenger()
{
  delete fUpdateCmd;
  delete fUseGeometryMacCmd;
  delete fManualFillterPosRatioZCmd;
  delete fManualFillterPosRatioYCmd;
  delete fManualFillterRatioZCmd;
  delete fManualFillterRatioYCmd;
  delete fManualCrystalSizeYCmd;
  delete fManualCrystalSizeCmd;
  delete fManualCrystalGapCmd;
  delete fManualArrayNzCmd;
  delete fManualArrayNyCmd;
  delete fManualArrayNxCmd;
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
  else if (command == fUseGeometryMacCmd)
    fDetector->SetUseGeometryMac(fUseGeometryMacCmd->GetNewBoolValue(newValue));
  else if (command == fManualArrayNxCmd)
    fDetector->SetManualArrayNx(fManualArrayNxCmd->GetNewIntValue(newValue));
  else if (command == fManualArrayNyCmd)
    fDetector->SetManualArrayNy(fManualArrayNyCmd->GetNewIntValue(newValue));
  else if (command == fManualArrayNzCmd)
    fDetector->SetManualArrayNz(fManualArrayNzCmd->GetNewIntValue(newValue));
  else if (command == fManualCrystalGapCmd)
    fDetector->SetManualCrystalGap(fManualCrystalGapCmd->GetNewDoubleValue(newValue));
  else if (command == fManualCrystalSizeCmd)
    fDetector->SetManualCrystalSize(fManualCrystalSizeCmd->GetNewDoubleValue(newValue));
  else if (command == fManualCrystalSizeYCmd)
    fDetector->SetManualCrystalSizeY(fManualCrystalSizeYCmd->GetNewDoubleValue(newValue));
  else if (command == fManualFillterRatioYCmd)
    fDetector->SetManualFillterRatioY(fManualFillterRatioYCmd->GetNewDoubleValue(newValue));
  else if (command == fManualFillterRatioZCmd)
    fDetector->SetManualFillterRatioZ(fManualFillterRatioZCmd->GetNewDoubleValue(newValue));
  else if (command == fManualFillterPosRatioYCmd)
    fDetector->SetManualFillterPosRatioY(fManualFillterPosRatioYCmd->GetNewDoubleValue(newValue));
  else if (command == fManualFillterPosRatioZCmd)
    fDetector->SetManualFillterPosRatioZ(fManualFillterPosRatioZCmd->GetNewDoubleValue(newValue));
  else if (command == fUpdateCmd) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
    G4cout << "Geometry updated with current parameters." << G4endl;
  }
}

}  // namespace B1

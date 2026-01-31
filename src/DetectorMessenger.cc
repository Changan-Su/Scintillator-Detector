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

  fUpdateCmd = new G4UIcmdWithoutParameter("/detector/update", this);
  fUpdateCmd->SetGuidance("Update geometry: rebuild detector with current parameters.");
  fUpdateCmd->SetGuidance("Use this after changing array size, gap, or crystal dimensions.");
}

DetectorMessenger::~DetectorMessenger()
{
  delete fUpdateCmd;
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
  else if (command == fUpdateCmd) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
    G4cout << "Geometry updated with current parameters." << G4endl;
  }
}

}  // namespace B1

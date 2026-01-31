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
/// \file DetectorMessenger.hh
/// \brief Definition of the B1::DetectorMessenger class

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithoutParameter;

namespace B1
{

class DetectorConstruction;

class DetectorMessenger : public G4UImessenger
{
  public:
    explicit DetectorMessenger(DetectorConstruction* det);
    ~DetectorMessenger() override;

    void SetNewValue(G4UIcommand* command, G4String newValue) override;

  private:
    DetectorConstruction* fDetector = nullptr;
    G4UIdirectory* fDetDir = nullptr;
    G4UIcmdWithAnInteger* fArrayNxCmd = nullptr;
    G4UIcmdWithAnInteger* fArrayNyCmd = nullptr;
    G4UIcmdWithAnInteger* fArrayNzCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fCrystalGapCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fCrystalSizeCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fCrystalSizeYCmd = nullptr;
    G4UIcmdWithADouble* fFillterRatioYCmd = nullptr;
    G4UIcmdWithADouble* fFillterRatioZCmd = nullptr;
    G4UIcmdWithADouble* fFillterPosRatioYCmd = nullptr;
    G4UIcmdWithADouble* fFillterPosRatioZCmd = nullptr;
    G4UIcmdWithoutParameter* fUpdateCmd = nullptr;
};

}  // namespace B1

#endif

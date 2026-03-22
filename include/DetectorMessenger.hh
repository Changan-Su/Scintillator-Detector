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
class G4UIcmdWithABool;
class G4UIcmdWithAString;

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
    G4UIcmdWithABool* fUseGeometryMacCmd = nullptr;
    G4UIcmdWithAnInteger* fManualArrayNxCmd = nullptr;
    G4UIcmdWithAnInteger* fManualArrayNyCmd = nullptr;
    G4UIcmdWithAnInteger* fManualArrayNzCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fManualCrystalGapCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fManualCrystalSizeCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fManualCrystalSizeYCmd = nullptr;
    G4UIcmdWithADouble* fManualFillterRatioYCmd = nullptr;
    G4UIcmdWithADouble* fManualFillterRatioZCmd = nullptr;
    G4UIcmdWithADouble* fManualFillterPosRatioYCmd = nullptr;
    G4UIcmdWithADouble* fManualFillterPosRatioZCmd = nullptr;
    G4UIcmdWithADouble* fSurfaceSigmaCmd = nullptr;
    G4UIdirectory* fResultsDir = nullptr;
    G4UIcmdWithAString* fResultsPrefixCmd = nullptr;
};

}  // namespace B1

#endif

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
/// \file B1/include/DetectorConstruction.hh
/// \brief Definition of the B1::DetectorConstruction class

#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class DetectorMessenger;

namespace B1
{

/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction() = default;
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;
    G4VPhysicalVolume* fPhysCrystal = nullptr;  // Physical volume for the crystal
    G4VPhysicalVolume* fLogicSiPM = nullptr;      // Physical volume for the PMT

    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
    G4LogicalVolume* GetlogicSiPM() const { return flogicSiPM; }  // Getter for PMT logical volume

    G4int GetCrystal_gap() const { return fCrystal_gap; } // Getter for crystal gap
    G4double GetCrystal_x() const { return fCrystal_x; }
    G4double GetCrystal_y() const { return fCrystal_y; }
    G4double GetCrystal_z() const { return fCrystal_z; } // Getter for
    G4int GetCrystal_nx() const { return fCrystal_nx; } // Getter for number of crystals in one dimension
    G4int GetCrystal_ny() const { return fCrystal_ny; } // Getter for number of crystals in the other dimension
    G4int GetCrystal_nz() const { return fCrystal_nz; } // Getter for number of crystals in height
    G4double Getcrystal_l() const { return fcrystal_l; } // Getter for crystal length
    G4double Getcrystal_ly() const { return fcrystal_ly; } // Getter for crystal width


      // —— 新增：供 Messenger 写入的 setter（做最基本的保护）
    void SetArrayNx(G4int v) { fPar_nx = std::max(1, v); }
    void SetArrayNy(G4int v) { fPar_ny = std::max(1, v); }
    void SetArrayNz(G4int v) { fPar_nz = std::max(1, v); }

  protected:
    G4LogicalVolume* fScoringVolume = nullptr;
    G4LogicalVolume* flogicSiPM = nullptr;  // PMT logical volume for optical photon tracking

  private:

    G4int fPar_nx = 11;
    G4int fPar_ny = 4;
    G4int fPar_nz = 4;

    G4int fCrystal_gap ; // Gap between crystals
    G4int fCrystal_nx; // Number of crystals in one dimension
    G4int fCrystal_ny; // Number of crystals in the other dimension
    G4int fCrystal_nz; // Number of crystals in height
    G4double fcrystal_l; // Length of the crystal
    G4double fcrystal_ly; // Width of the crystal
    G4double fCrystal_x; // Crystal size in x direction
    G4double fCrystal_y; // Crystal size in y direction
    G4double fCrystal_z; // Crystal size in z direction

    // Messenger 指针
    DetectorMessenger* fMessenger = nullptr;
  };

}  // namespace B1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

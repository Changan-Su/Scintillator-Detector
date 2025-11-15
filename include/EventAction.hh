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
/// \file B1/include/EventAction.hh
/// \brief Definition of the B1::EventAction class

#ifndef B1EventAction_h
#define B1EventAction_h 1
// #pragma message(">>> CORRECT EventAction.hh is: " __FILE__)

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "HistoManager.hh"



class G4Event;

namespace B1
{

class RunAction;

/// Event action class

class EventAction : public G4UserEventAction
{
  public:
  EventAction(RunAction* runAction, HistoManager* histo); //  改为包含两个参数
  ~EventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;
    // 计数 optical photon
    


    void AddEdep(G4double edep) { fEdep += edep; }
    // void AddPhotonDepth(G4int depth) {fHistoManager->FillPhotonDepth(depth);} 
    std::vector<G4int> fPhotonDOIs;  // 所有光子对应的 depth
    void AddPhotonDepth(G4int depth) { fPhotonDOIs.push_back(depth); }  

    void AddAbsorption(G4double edep, G4double trackL) {
      fEnergyAbs += edep;
      fTrackLAbs += trackL;
      }
    void AddPhotonInDetector(G4int count) { fPhotonCountDetector += count ;} // 统计进入探测器的光子数量
    void AddPhotonInPMT(G4int count) { fPhotonCountPMT += count ;} // 统计进入 PMT 的光子数量
    
    // void AddPhotonFuck(G4int count) { fPhotonCountLeft += count ;} // 统计生成的光子数量
    G4int FuckLeft = 1;
    void AddPhotonLeft(G4int count) { fPhotonCountLeft += count; } // 统计进入左侧探测器的光子数量
    void AddPhotonRight(G4int count) { fPhotonCountRight += count; }

    G4int GetPhotonCountDetector() const { return fPhotonCountDetector; }// 获取进入探测器的光子数量
    G4int GetPhotonCountPMT() const { return fPhotonCountPMT; }// 获取进入 PMT 的光子数量
    G4int GetPhotonCountLeft() const { return fPhotonCountLeft; }// 获取进入左侧探测器的光子数量
    G4int GetPhotonCountRight() const { return fPhotonCountRight; }// 获取进入右侧探测器的光子数量
    void AddPhotonGenerated(G4int n = 1) { fPhotonCountGenerated += n; }// 统计生成的光子数量
    G4int GetPhotonGenerated() const { return fPhotonCountGenerated; }// 获取生成的光子数量
    G4int Fuckrigt = 2;
    void AddPhotonLeftAt(G4int iy, G4int iz)  { fLeftPerRod[RodIndex(iy,iz)]++; }
    void AddPhotonRightAt(G4int iy, G4int iz) { fRightPerRod[RodIndex(iy,iz)]++; }
    

    private:
    RunAction* fRunAction = nullptr;
    G4double fEdep = 0.;
    G4double fEnergyAbs = 0.;
    G4double fTrackLAbs = 0.;
    G4int fPrintModulo = 1000;
    G4bool fApplySmearing = false;
    G4int fPhotonCountGenerated; 
    G4int fPhotonCountDetector = 0; // 统计进入探测器的光子数量
    G4int fPhotonCountPMT = 0; // 统计进入 PMT 的光子数量
    G4int fPhotonFuck = 0;
    G4int fPhotonCountLeft = 0;
    G4int fPhotonCountRight = 0;


    std::vector<int> fLeftPerRod, fRightPerRod;
    G4int fNy=0, fNz=0;
    inline int RodIndex(int iy, int iz) const { return iz*fNy + iy; }


    

    HistoManager* fHistoManager = nullptr;
};




}  // namespace B1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

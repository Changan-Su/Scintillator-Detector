#include "MyPhysicsList.hh"
#include "G4DecayPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalParameters.hh" 
#include "G4Scintillation.hh"
#include "G4ProcessManager.hh"
#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4Scintillation.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"

MyPhysicsList::MyPhysicsList()
{
    SetVerboseLevel(1);

    AddTransportation();
    RegisterPhysics(new G4HadronPhysicsFTFP_BERT());
    RegisterPhysics(new G4DecayPhysics());
    RegisterPhysics(new G4EmLivermorePhysics());
    // ReplacePhysics(new G4EmStandardPhysics_option4());

    auto opticalPhysics = new G4OpticalPhysics();
    RegisterPhysics(opticalPhysics);
    auto params = G4OpticalParameters::Instance();
    params->SetWLSTimeProfile("delta");

    params->SetScintTrackSecondariesFirst(true);
    params->SetCerenkovTrackSecondariesFirst(true); // 开启 Cerenkov 光子跟踪
    params->SetCerenkovMaxPhotonsPerStep(0);
    params->SetCerenkovMaxBetaChange(10.0);
    params->SetScintByParticleType(true);          // 按粒子种类区分闪烁光
    params->SetScintStackPhotons(true);             // Scintillation光子叠加到堆栈


    auto atomDeexcitation = new G4UAtomicDeexcitation();
    atomDeexcitation->SetFluo(true);//原子退激发
    atomDeexcitation->SetAuger(true);
    atomDeexcitation->SetPIXE(true);
    
    G4int FuckU = 1;
    G4LossTableManager::Instance()->SetAtomDeexcitation(atomDeexcitation);
    //Fuck you
}




void MyPhysicsList::SetCuts()
{
    G4VUserPhysicsList::SetCuts();
}

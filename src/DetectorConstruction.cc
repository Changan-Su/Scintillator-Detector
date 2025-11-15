
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trd.hh" 

#include "G4Region.hh"
#include "G4ProductionCuts.hh"

// #include "DetectorMessenger.hh"  

namespace B1
{

  
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  G4VPhysicalVolume* DetectorConstruction::Construct()
  {

    G4double Surface_Sigma = 0.5;

    // Get nist material manager
    G4NistManager* nist = G4NistManager::Instance();
  
    // Envelope parameters
    //
    G4double env_sizeXY = 50 * cm, env_sizeZ = 50 * cm;
    G4Material* env_mat = nist->FindOrBuildMaterial("G4_Galactic");
    auto env_mt = new G4MaterialPropertiesTable();
    std::vector<G4double> energy_env = {2.0 * eV, 3.5 * eV};
    std::vector<G4double> rindex_env = {1.0, 1.0};
       // std::vector<G4double> abslen_env = {1000 * cm, 1000 * cm}; 
    env_mt->AddProperty("RINDEX", energy_env, rindex_env);
    // env_mt->AddProperty("ABSLENGTH", energy_env, abslen_env);
    env_mat->SetMaterialPropertiesTable(env_mt);

    G4Material* env_mat_air = nist -> FindOrBuildMaterial("G4_AIR");
    auto env_mt_air = new G4MaterialPropertiesTable();
    std::vector<G4double> energy_env_air = {2.0 * eV, 3.5 * eV};
    std::vector<G4double> rindex_env_air = {1.0003, 1.0003};
    env_mt_air->AddProperty("RINDEX", energy_env_air, rindex_env_air);
    env_mat_air->SetMaterialPropertiesTable(env_mt_air);


 
 
    // Option to switch on/off checking of volumes overlaps
    //
    G4bool checkOverlaps = true;
  
    //
    // World
    //
    
    G4double world_sizeXY = env_sizeXY;
    G4double world_sizeZ =  env_sizeZ;
    G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
   
    
  
    auto solidWorld =
      new G4Box("World",  // its name
                0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);  // its size
  
    auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
                                          world_mat,  // its material
                                          "World");  // its name
  
    auto physWorld = new G4PVPlacement(nullptr,  // no rotation
                                       G4ThreeVector(),  // at (0,0,0)
                                       logicWorld,  // its logical volume
                                       "World",  // its name
                                       nullptr,  // its mother  volume
                                       false,  // no boolean operation
                                       0,  // copy number
                                       checkOverlaps);  // overlaps checking
  
    //
    // Envelope
    //
    auto solidEnv = new G4Box("Envelope",  // its name
                              0.5 * env_sizeXY, 0.5 * env_sizeXY, 0.5 * env_sizeZ);  // its size
  
    auto logicEnv = new G4LogicalVolume(solidEnv,  // its solid
                                        env_mat_air,  // its material
                                        "Envelope");  // its name
  
    auto physEnv = new G4PVPlacement(nullptr,  // no rotation
                      G4ThreeVector(),  // at (0,0,0)
                      logicEnv,  // its logical volume
                      "Envelope",  // its name
                      logicWorld,  // its mother  volume
                      false,  // no boolean operation
                      0,  // copy number
                      checkOverlaps);  // overlaps checking
//=========================================================================================

    // G4Element* elNa = G4NistManager::Instance()->FindOrBuildElement("Na");
    // G4Element* elI  = G4NistManager::Instance()->FindOrBuildElement("I");
    // G4double density = 3.67 * g/cm3;
    // G4Material* NaI_Tl = new G4Material("NaI_Tl", density, 2);
    // NaI_Tl->AddElement(elNa, 1);
    // NaI_Tl->AddElement(elI, 1);

    // // 材料属性表用LXe风格
    // std::vector<G4double> nai_Energy = {2.07 * eV, 2.34 * eV, 2.62 * eV, 2.89 * eV, 3.10 * eV};
    // std::vector<G4double> nai_SCINT = {1.0, 1.0, 1.0, 1.0, 1.0}; // 保持能量范围覆盖NaI发射区间
    // std::vector<G4double> nai_RIND = {1.85, 1.85, 1.85, 1.85, 1.85};
    // // std::vector<G4double> nai_ABSL = {38. * cm, 38. * cm, 38. * cm, 38. * cm, 38. * cm};
    // std::vector<G4double> nai_ABSL = {40. * cm, 40. * cm, 40. * cm, 40. * cm, 40. * cm}; // 吸收长度设置为100cm

    // auto nai_mt = new G4MaterialPropertiesTable();
    // nai_mt->AddProperty("SCINTILLATIONCOMPONENT1", nai_Energy, nai_SCINT);
    // nai_mt->AddProperty("SCINTILLATIONCOMPONENT2", nai_Energy, nai_SCINT); // 没slow就全0也行
    // nai_mt->AddProperty("RINDEX", nai_Energy, nai_RIND);
    // nai_mt->AddProperty("ABSLENGTH", nai_Energy, nai_ABSL);
    // nai_mt->AddConstProperty("SCINTILLATIONYIELD", 38000. / MeV); // 或者先试12000看看能否产光
    // nai_mt->AddConstProperty("RESOLUTIONSCALE", 1.0);
    // nai_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 250. * ns);
    // nai_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 0. * ns);
    // nai_mt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
    // nai_mt->AddConstProperty("SCINTILLATIONYIELD2", 0.0);
    // NaI_Tl->SetMaterialPropertiesTable(nai_mt);


    //======================= YSO Material Definition =======================
    G4Element* elY = G4NistManager::Instance()->FindOrBuildElement("Y");
    G4Element* elSi = G4NistManager::Instance()->FindOrBuildElement("Si");
    G4Element* elO = G4NistManager::Instance()->FindOrBuildElement("O");
    G4double density_YSO = 4.45 * g/cm3;
    G4Material* YSO = new G4Material("YSO", density_YSO, 3);
    YSO->AddElement(elY, 2);
    YSO->AddElement(elSi, 1);
    YSO->AddElement(elO, 5);

    // 发射峰在420 nm => 2.95 eV 附近
    std::vector<G4double> yso_Energy = {2.07 * eV, 2.34 * eV, 2.62 * eV, 2.89 * eV, 3.10 * eV};
    std::vector<G4double> yso_SCINT = {1.0, 1.0, 1.0, 1.0, 1.0}; // 简单设定恒定发光强度
    G4int Fukkkk = 0.6 *cm;
    std::vector<G4double> yso_RIND = {1.8, 1.8, 1.8, 1.8, 1.8};  // 折射率 ~1.8
    std::vector<G4double> yso_ABSL = {40. * cm, 40. * cm, 40. * cm, 40. * cm, 40. * cm}; // 吸收长度

    auto yso_mt = new G4MaterialPropertiesTable();
    yso_mt->AddProperty("SCINTILLATIONCOMPONENT1", yso_Energy, yso_SCINT);
    yso_mt->AddProperty("SCINTILLATIONCOMPONENT2", yso_Energy, yso_SCINT); // 没 slow 成分可设为 0

    yso_mt->AddProperty("RINDEX", yso_Energy, yso_RIND);
    yso_mt->AddProperty("ABSLENGTH", yso_Energy, yso_ABSL);
    yso_mt->AddConstProperty("SCINTILLATIONYIELD", 24000. / MeV); // 光产额
    G4int Fuck666 = 1*cm;
    yso_mt->AddConstProperty("RESOLUTIONSCALE", 1.0);
    yso_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 56. * ns); // 衰减时间
    yso_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 0. * ns);  // 无慢分量
    yso_mt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
    yso_mt->AddConstProperty("SCINTILLATIONYIELD2", 0.0);
    YSO->SetMaterialPropertiesTable(yso_mt);

    //LocCry
    
      G4double Crystal_gap = 0.1 * mm;
      G4int Crystal_nx = 11;
      G4int Crystal_ny = 1;
      G4int Crystal_nz = 1;
      G4double crystal_l = 3 * mm; // 单个晶体的边长
      // G4double FuckGeant4 = 6 *cm;//是的，你没看错，加了这一行代码就不报错了
      G4double crystal_ly = 1 * crystal_l;

      G4double Crystal_x = crystal_l * Crystal_nx + Crystal_gap * (Crystal_nx - 1);
      G4double Crystal_y = crystal_ly * Crystal_ny + Crystal_gap * (Crystal_ny - 1);
      G4double Crystal_z = crystal_l * Crystal_nz + Crystal_gap * (Crystal_nz - 1);
      
      auto solidCrystal = new G4Box("Crystal", crystal_l / 2, crystal_ly / 2, crystal_l / 2);
      auto logicCrystal = new G4LogicalVolume(solidCrystal, YSO, "Crystal");


      //Crystal Fillter
      G4double Fillter_x = Crystal_gap;
      G4double Fillter_y = crystal_ly * 0.3;
      G4double Fillter_z = crystal_l * 0.2;
      G4bool IfFillter = true;

      auto solidFillter = new G4Box("Fillter", Fillter_x/2, Fillter_y/2, Fillter_z/2);
      auto logicFillter = new G4LogicalVolume(solidFillter, YSO, "Fillter");


      //SiPM Mat
      G4Material* SiPM_mat = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
      G4double sipm_l = crystal_l ; // SiPM size
      G4double sipm_l_ratio = 0.1;
      G4int SiPm_nz = Crystal_z / sipm_l ;//注意Gap
      G4int SiPm_nx = Crystal_x / sipm_l ;
      G4Box* solidSiPM = new G4Box("SiPM", sipm_l * sipm_l_ratio/ 2, sipm_l / 2, sipm_l / 2);
 
      
      std::vector<G4double> energy = {2.0*eV, 3.5*eV};
      std::vector<G4double> rindex_sipm = {1.5, 1.5};  // 硅的折射率 ~1.5
      G4LogicalVolume* logicSiPM = new G4LogicalVolume(solidSiPM, SiPM_mat, "SiPM");
      auto sipm_mt = new G4MaterialPropertiesTable();
      sipm_mt->AddProperty("RINDEX", energy, rindex_sipm);
      SiPM_mat->SetMaterialPropertiesTable(sipm_mt);


      for (int iz = 0; iz < Crystal_nz; ++iz)
      {
      G4double posZ = -Crystal_z/2 + iz * (crystal_l + Crystal_gap) + crystal_l / 2;
        for (int iy = 0; iy < Crystal_ny; ++iy)
        {
          G4double posY = -Crystal_y/2 + iy * (crystal_ly + Crystal_gap) + crystal_ly / 2;
          for (int ix = 0; ix < Crystal_nx; ++ix)
          {
            G4double posX = -Crystal_x/2 + ix * (crystal_l + Crystal_gap) + crystal_l / 2;
            G4ThreeVector pos_crystal = G4ThreeVector(posX, posY, posZ);
            G4int copyNo = 0;
            auto physCrystal = new G4PVPlacement(nullptr,  // no rotation
                                                  pos_crystal,  // position
                                                  logicCrystal,  // its logical volume
                                                  "Unit_crystal",  // its name
                                                  logicEnv,  // its mother volume
                                                  false,  // no boolean operation
                                                  copyNo,  // copy number
                                                  checkOverlaps);  // overlaps checking
            if(ix != 0 && IfFillter)
            {
              G4double pos_fillter_x = -Crystal_x/2 + (crystal_l + Fillter_x) * ix - Fillter_x/2;
              G4ThreeVector pos_fillter = G4ThreeVector(pos_fillter_x, posY, posZ);
              auto physFillter = new G4PVPlacement(nullptr,
                                                   pos_fillter,
                                                   logicFillter,
                                                   "Unit_fillter",
                                                   logicEnv,
                                                   false,
                                                   copyNo,
                                                   checkOverlaps);
            }
            // auto CrystalToAir = new G4OpticalSurface("CrystalToAir");
            // CrystalToAir->SetType(dielectric_dielectric);
            // CrystalToAir->SetModel(unified);
            // CrystalToAir->SetFinish(polished);          // 粗糙侧面
            // CrystalToAir->SetSigmaAlpha(Surface_Sigma);

            // std::string surfName = "Crystal-Air-" +
            //                       std::to_string(ix) + "-" +
            //                       std::to_string(iy) + "-" +
            //                       std::to_string(iz);

            // new G4LogicalBorderSurface(surfName,
            //                           physCrystal, physEnv, CrystalToAir);

          }

      G4double SiPM_PosX_l = -Crystal_x/2 - sipm_l * sipm_l_ratio/2;
      G4double SiPM_PosX_r = Crystal_x/2 + sipm_l * sipm_l_ratio/2;
      //FuckSiPM
      // G4double SiPM_PosX_l = -Crystal_x/2 - sipm_l * sipm_l_ratio/2 - 0.1 * mm;
      // G4double SiPM_PosX_r = Crystal_x/2 + sipm_l * sipm_l_ratio/2 + 0.1 * mm;
      G4cout <<"DEBUG:SIPM_POSR"<<SiPM_PosX_r<<G4endl;
      G4cout <<"DEBUG:CRYSTAL_X"<<Crystal_x/2<<G4endl;

      G4double SiPM_PosY = posY;
      G4double SiPM_PosZ = posZ;
      G4ThreeVector SiPM_Pos_l = G4ThreeVector(SiPM_PosX_l, SiPM_PosY, SiPM_PosZ);
      G4ThreeVector SiPM_Pos_r = G4ThreeVector(SiPM_PosX_r, SiPM_PosY, SiPM_PosZ);
      new G4PVPlacement(
              nullptr,
              SiPM_Pos_l,
              logicSiPM,
              "SiPM_Left",
              logicEnv,
              false,
              iz * 1000000 + iy * 100 + 1,
              checkOverlaps
            );
      new G4PVPlacement(
              nullptr,
              SiPM_Pos_r,
              logicSiPM,
              "SiPM_Right",
              logicEnv,
              false,
              (iz) * 1000000 + iy * 100 + 2,
              checkOverlaps
            );


        }
      }



      fScoringVolume = logicCrystal;
      fCrystal_gap = Crystal_gap; // Store crystal gap
      fCrystal_nx = Crystal_nx; // Store number of crystals in one dimension
      fCrystal_ny = Crystal_ny; // Store number of crystals in the other dimension
      fCrystal_nz = Crystal_nz; // Store number of crystals in height
      fcrystal_l = crystal_l; // Store crystal length
      fcrystal_ly = crystal_ly; // Store crystal width
      fCrystal_x = Crystal_x; // Store crystal size in x direction
      fCrystal_y = Crystal_y; // Store crystal size in y direction
      fCrystal_z = Crystal_z; // Store crystal size in z direction


      // //Crystal Optical Surface
      // G4OpticalSurface* crystalsurface = new G4OpticalSurface("CrystalSurface");
      // crystalsurface->SetType(dielectric_dielectric);
      // crystalsurface->SetModel(unified);
      // crystalsurface->SetFinish(polished);
      // crystalsurface->SetSigmaAlpha(Surface_Sigma);
      // G4cout << "DEBUG: CrystalSurface SigmaAlpha = "
      //  << crystalsurface->GetSigmaAlpha() << G4endl;
      // new G4LogicalSkinSurface("CrystalSurface",logicCrystal, crystalsurface);






      //SurfSiP



       // 1) 皮肤光学表面，只创建一次

      // G4OpticalSurface* SiPM_Surf = new G4OpticalSurface("SiPMSkinSurface");
      // SiPM_Surf->SetType(dielectric_dielectric);
      // SiPM_Surf->SetModel(unified);
      // SiPM_Surf->SetFinish(polished);

      // // 2) 贴到整个 logicSiPM 上
      // new G4LogicalSkinSurface("SiPMSkinSurface",logicSiPM,SiPM_Surf);




      // //Optical Surface between air and crystal
      // G4OpticalSurface* Air_Crystal_Surf = new G4OpticalSurface("Air_Crystal_Surface");
      // Air_Crystal_Surf -> SetType(dielectric_dielectric);
      // Air_Crystal_Surf -> SetModel(unified);
      // Air_Crystal_Surf->SetFinish(polished);       // 或 polished
      // // Air_Crystal_Surf->SetSigmaAlpha(0);      // 需要粗糙度才设
      // // new G4LogicalBorderSurface("Crystal-Air", physCrystal, physEnv, Air_Crystal_Surf);
      // new G4LogicalSkinSurface("Air_Crystal_Surface", logicCrystal, Air_Crystal_Surf);

      //3 Double-end SiPM
 
      
      
      // // Position SiPM on the crystal
      // for (int iz = 0;iz < SiPm_nz ; ++iz)
      //   {
      //     for (int ix = 0 ;ix < SiPm_nx; ++ix)
      //     {
      //       G4double PosX = -Crystal_x/2 + ix * sipm_l;
      //       G4double PosY = Crystal_y/2 + sipm_l/2 ;
      //       G4double PosZ = -Crystal_z/2 + iz * sipm_l;

      //       G4ThreeVector SiPm_Pos_Top = G4ThreeVector(PosX,PosY,PosZ);
      //       G4ThreeVector SiPm_Pos_Bottom = G4ThreeVector(PosX,-PosY,PosZ);

      //       new G4PVPlacement(
      //         nullptr,
      //         SiPm_Pos_Top,
      //         logicSiPM,
      //         "SiPM_Top",
      //         logicEnv,
      //         false,
      //         iz*100000+ix*10+1
      //       );
      //       new G4PVPlacement(
      //         nullptr,
      //         SiPm_Pos_Bottom,
      //         logicSiPM,
      //         "SiPM_Bottom",
      //         logicEnv,
      //         false,
      //         iz*100000+ix*10+2
      //       );
      //     }
      //   }
      flogicSiPM = logicSiPM;
    return physWorld;
  }
}

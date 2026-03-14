
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
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
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trd.hh" 

#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include <string>
#include <cstdio>

namespace B1
{

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  DetectorConstruction::DetectorConstruction()
  {
    fMessenger = new DetectorMessenger(this);
  }

  DetectorConstruction::~DetectorConstruction()
  {
    delete fMessenger;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4VPhysicalVolume* DetectorConstruction::Construct()
  {
    G4double Surface_Sigma = fSurfaceSigma;

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


    // //======================= YSO Material Definition =======================
    // G4Element* elY = G4NistManager::Instance()->FindOrBuildElement("Y");
    // G4Element* elSi = G4NistManager::Instance()->FindOrBuildElement("Si");
    // G4Element* elO = G4NistManager::Instance()->FindOrBuildElement("O");
    // G4double density_YSO = 4.45 * g/cm3;
    // G4Material* YSO = new G4Material("YSO", density_YSO, 3);
    // YSO->AddElement(elY, 2);
    // YSO->AddElement(elSi, 1);
    // YSO->AddElement(elO, 5);

    // // 发射峰在420 nm => 2.95 eV 附近
    // std::vector<G4double> yso_Energy = {2.07 * eV, 2.34 * eV, 2.62 * eV, 2.89 * eV, 3.10 * eV};
    // std::vector<G4double> yso_SCINT = {1.0, 1.0, 1.0, 1.0, 1.0}; // 简单设定恒定发光强度
    // G4int Fukkkk = 0.6 *cm;
    // std::vector<G4double> yso_RIND = {1.8, 1.8, 1.8, 1.8, 1.8};  // 折射率 ~1.8
    // std::vector<G4double> yso_ABSL = {1.5 * cm, 1.5 * cm, 1.5 * cm, 1.5 * cm, 1.5 * cm}; // 吸收长度

    // auto yso_mt = new G4MaterialPropertiesTable();
    // yso_mt->AddProperty("SCINTILLATIONCOMPONENT1", yso_Energy, yso_SCINT);
    // yso_mt->AddProperty("SCINTILLATIONCOMPONENT2", yso_Energy, yso_SCINT); // 没 slow 成分可设为 0

    // yso_mt->AddProperty("RINDEX", yso_Energy, yso_RIND);
    // yso_mt->AddProperty("ABSLENGTH", yso_Energy, yso_ABSL);
    // yso_mt->AddConstProperty("SCINTILLATIONYIELD", 24000. / MeV); // 光产额
    // G4int Fuck666 = 1*cm;
    // yso_mt->AddConstProperty("RESOLUTIONSCALE", 1.0);
    // yso_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 56. * ns); // 衰减时间
    // yso_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 0. * ns);  // 无慢分量
    // yso_mt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
    // yso_mt->AddConstProperty("SCINTILLATIONYIELD2", 0.0);
    // YSO->SetMaterialPropertiesTable(yso_mt);


    // GAGG scintillator: Gd3Al2Ga3O12
    G4Element* elGd = G4NistManager::Instance()->FindOrBuildElement("Gd");
    G4Element* elAl = G4NistManager::Instance()->FindOrBuildElement("Al");
    G4Element* elGa = G4NistManager::Instance()->FindOrBuildElement("Ga");
    G4Element* elO = G4NistManager::Instance()->FindOrBuildElement("O");
    G4double density_GAGG = 6.6 * g/cm3;
    G4Material* GAGG = new G4Material("GAGG", density_GAGG, 4);
    GAGG->AddElement(elGd, 3);
    GAGG->AddElement(elAl, 2);
    GAGG->AddElement(elGa, 3);
    GAGG->AddElement(elO, 12);

    // GAGG-HL optical/scint properties: n=1.91, 54k ph/MeV, 150 ns decay, ~530 nm peak
    std::vector<G4double> gagg_Energy = {2.07 * eV, 2.34 * eV, 2.62 * eV, 2.89 * eV, 3.10 * eV};
    std::vector<G4double> gagg_SCINT = {1.0, 1.0, 1.0, 1.0, 1.0};
    std::vector<G4double> gagg_RIND = {1.91, 1.91, 1.91, 1.91, 1.91};
    std::vector<G4double> gagg_ABSL = {2.5 * cm, 2.5 * cm, 2.5 * cm, 2.5 * cm, 2.5 * cm};
    auto gagg_mt = new G4MaterialPropertiesTable();
    gagg_mt->AddProperty("SCINTILLATIONCOMPONENT1", gagg_Energy, gagg_SCINT);
    gagg_mt->AddProperty("SCINTILLATIONCOMPONENT2", gagg_Energy, gagg_SCINT);
    gagg_mt->AddProperty("RINDEX", gagg_Energy, gagg_RIND);
    gagg_mt->AddProperty("ABSLENGTH", gagg_Energy, gagg_ABSL);
    gagg_mt->AddConstProperty("SCINTILLATIONYIELD", 54000. / MeV);
    gagg_mt->AddConstProperty("RESOLUTIONSCALE", 1.0);
    gagg_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 150. * ns);
    gagg_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 0. * ns);
    gagg_mt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
    gagg_mt->AddConstProperty("SCINTILLATIONYIELD2", 0.0);
    GAGG->SetMaterialPropertiesTable(gagg_mt);

    // Crystal array params: true = from geometry.mac/UI, false = manual below
    G4int Crystal_nx, Crystal_ny, Crystal_nz;
    G4double Crystal_gap, crystal_l, crystal_ly;
    G4double Fillter_Gap_Ratio_Y, Fillter_Gap_Ratio_Z, Fillter_Gap_PosRatio_Y, Fillter_Gap_PosRatio_Z;
    G4bool fUseGeometryMac = false;
    if (fUseGeometryMac) {
      Crystal_nx = fPar_nx;
      Crystal_ny = fPar_ny;
      Crystal_nz = fPar_nz;
      Crystal_gap = fCrystal_gap * mm;
      crystal_l = fcrystal_l * mm;
      crystal_ly = fcrystal_ly * mm;
      Fillter_Gap_Ratio_Y = fFillter_Gap_Ratio_Y;
      Fillter_Gap_Ratio_Z = fFillter_Gap_Ratio_Z;
      Fillter_Gap_PosRatio_Y = fFillter_Gap_PosRatio_Y;
      Fillter_Gap_PosRatio_Z = fFillter_Gap_PosRatio_Z;
    } else {
      Crystal_nx = 1;
      Crystal_ny = 1;
      Crystal_nz = 1;
      Crystal_gap = 0.0 * mm;
      crystal_l = 25 * mm;
      crystal_ly = 25 * mm;
      Fillter_Gap_Ratio_Y = 0.3;
      Fillter_Gap_Ratio_Z = 1.0;
      Fillter_Gap_PosRatio_Y = 0.7;
      Fillter_Gap_PosRatio_Z = 0.0;
    }

    // Array total size; single crystal solid and logical
    G4double Crystal_x = crystal_l * Crystal_nx + Crystal_gap * (Crystal_nx - 1);
    G4double Crystal_y = crystal_ly * Crystal_ny + Crystal_gap * (Crystal_ny - 1);
    G4double Crystal_z = crystal_l * Crystal_nz + Crystal_gap * (Crystal_nz - 1);
    G4Box* solidCrystal = new G4Box("Crystal", crystal_l / 2, crystal_ly / 2, crystal_l / 2);
    G4LogicalVolume* logicCrystal = new G4LogicalVolume(solidCrystal, GAGG, "Crystal");

    // Fillter (gap filler), only when Crystal_gap > 0
    G4double Fillter_x = Crystal_gap;
    G4double Fillter_y = crystal_ly * Fillter_Gap_Ratio_Y;
    G4double Fillter_z = crystal_l * Fillter_Gap_Ratio_Z;
    G4bool IfFillter = (Crystal_gap > 0.);
    G4LogicalVolume* logicFillter = nullptr;
    if (IfFillter) {
      auto solidFillter = new G4Box("Fillter", Fillter_x/2, Fillter_y/2, Fillter_z/2);
      logicFillter = new G4LogicalVolume(solidFillter, GAGG, "Fillter");
    }

    // SiPM: 6mm pixel, 0.2mm gap, 4x4 per face, 6 faces per crystal
    G4Material* SiPM_mat = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    G4double sipm_l = 6. * mm;
    G4double sipm_t = 0.6 * mm;
    G4double SiPm_gap = 0. * mm;
    const G4int SiPm_np = 4;
    G4Box* solidSiPM = new G4Box("SiPM", sipm_t / 2, sipm_l / 2, sipm_l / 2);
    std::vector<G4double> energy = {2.0*eV, 3.5*eV};
    std::vector<G4double> rindex_sipm = {1.5, 1.5};
    G4LogicalVolume* logicSiPM = new G4LogicalVolume(solidSiPM, SiPM_mat, "SiPM");
    auto sipm_mt = new G4MaterialPropertiesTable();
    sipm_mt->AddProperty("RINDEX", energy, rindex_sipm);
    SiPM_mat->SetMaterialPropertiesTable(sipm_mt);

    // Optical grease (硅脂, PDMS-based): n≈1.46, fills crystal-SiPM interface.
    // No explicit G4OpticalSurface is needed: Geant4 G4OpBoundaryProcess performs
    // automatic Fresnel reflection/refraction at any boundary where both materials
    // carry RINDEX (crystal n=1.91 → grease n=1.46 → SiPM n=1.5).
    G4Element* elC   = G4NistManager::Instance()->FindOrBuildElement("C");
    G4Element* elH   = G4NistManager::Instance()->FindOrBuildElement("H");
    G4Element* elSi  = G4NistManager::Instance()->FindOrBuildElement("Si");
    G4double   grease_t   = 0.1 * mm;
    G4Material* grease_mat = new G4Material("OpticalGrease", 1.06 * g/cm3, 4);
    grease_mat->AddElement(elC,  2);
    grease_mat->AddElement(elH,  6);
    grease_mat->AddElement(elSi, 1);
    grease_mat->AddElement(elO,  1);
    std::vector<G4double> grease_Energy = {2.0 * eV, 3.5 * eV};
    std::vector<G4double> grease_RIND   = {1.46, 1.46};
    std::vector<G4double> grease_ABSL   = {1000. * cm, 1000. * cm};
    auto grease_mt = new G4MaterialPropertiesTable();
    grease_mt->AddProperty("RINDEX",    grease_Energy, grease_RIND);
    grease_mt->AddProperty("ABSLENGTH", grease_Energy, grease_ABSL);
    grease_mat->SetMaterialPropertiesTable(grease_mt);
    // One grease slab per crystal face; two shapes because crystal_ly may differ from crystal_l.
    // Rotation convention mirrors SiPM: X-axis of the solid becomes the face-normal in world.
    //   Faces 0,1 (X±) and 4,5 (Z±): local YZ covers crystal_ly × crystal_l
    //   Faces 2,3 (Y±):               local YZ covers crystal_l  × crystal_l
    G4Box* solidGrease_XZ = new G4Box("Grease_XZ", grease_t / 2, crystal_ly / 2, crystal_l / 2);
    G4Box* solidGrease_Y  = new G4Box("Grease_Y",  grease_t / 2, crystal_l  / 2, crystal_l / 2);
    G4LogicalVolume* logicGrease_XZ = new G4LogicalVolume(solidGrease_XZ, grease_mat, "Grease_XZ");
    G4LogicalVolume* logicGrease_Y  = new G4LogicalVolume(solidGrease_Y,  grease_mat, "Grease_Y");

    // Keep rotation matrices alive during all placements.
    G4RotationMatrix* rotSiPM_mX = new G4RotationMatrix();
    rotSiPM_mX->rotateY(180. * deg);
    G4RotationMatrix* rotSiPM_pY = new G4RotationMatrix();
    rotSiPM_pY->rotateZ(90. * deg);
    G4RotationMatrix* rotSiPM_mY = new G4RotationMatrix();
    rotSiPM_mY->rotateZ(-90. * deg);
    G4RotationMatrix* rotSiPM_pZ = new G4RotationMatrix();
    rotSiPM_pZ->rotateY(-90. * deg);
    G4RotationMatrix* rotSiPM_mZ = new G4RotationMatrix();
    rotSiPM_mZ->rotateY(90. * deg);

    // Loop iz -> iy -> ix: place crystal, fillter, then 6 faces x 4x4 SiPM per crystal
    for (int iz = 0; iz < Crystal_nz; ++iz) {
      G4double posZ = -Crystal_z/2 + iz * (crystal_l + Crystal_gap) + crystal_l / 2;
      G4double Fillter_Pos_Z = posZ + crystal_l * Fillter_Gap_PosRatio_Z/2;
      for (int iy = 0; iy < Crystal_ny; ++iy) {
        G4double posY = -Crystal_y/2 + iy * (crystal_ly + Crystal_gap) + crystal_ly / 2;
        G4double Fillter_Pos_Y = posY + crystal_ly*Fillter_Gap_PosRatio_Y/2;
        for (int ix = 0; ix < Crystal_nx; ++ix) {
          G4double posX = -Crystal_x/2 + ix * (crystal_l + Crystal_gap) + crystal_l / 2;
          G4ThreeVector pos_crystal = G4ThreeVector(posX, posY, posZ);
          G4int copyNo = 0;
          new G4PVPlacement(nullptr, pos_crystal, logicCrystal, "Unit_crystal", logicEnv, false, copyNo, checkOverlaps);
          if (ix != 0 && IfFillter && logicFillter != nullptr) {
            G4double Fillter_Pos_X = -Crystal_x/2 + (crystal_l + Fillter_x) * ix - Fillter_x/2;
            G4ThreeVector pos_fillter = G4ThreeVector(Fillter_Pos_X, Fillter_Pos_Y, Fillter_Pos_Z);
            new G4PVPlacement(nullptr, pos_fillter, logicFillter, "Unit_fillter", logicEnv, false, copyNo, checkOverlaps);
          }

          // 6 faces of this crystal: 4x4 SiPM each, pitch 6.2 mm
          G4double pitch = sipm_l + SiPm_gap;
          G4int crystalId = ix * Crystal_ny * Crystal_nz + iy * Crystal_nz + iz;
          for (G4int face = 0; face < 6; ++face) {
            G4ThreeVector faceCenter;
            G4ThreeVector greaseCenter;
            G4RotationMatrix* rot = nullptr;
            G4LogicalVolume* logicGreaseFace = nullptr;
            if (face == 0) {//X+
              greaseCenter    = G4ThreeVector(posX + crystal_l/2 + grease_t/2, posY, posZ);
              faceCenter      = G4ThreeVector(posX + crystal_l/2 + grease_t + sipm_t/2, posY, posZ);
              logicGreaseFace = logicGrease_XZ;
            } else if (face == 1) {//X-
              greaseCenter    = G4ThreeVector(posX - crystal_l/2 - grease_t/2, posY, posZ);
              faceCenter      = G4ThreeVector(posX - crystal_l/2 - grease_t - sipm_t/2, posY, posZ);
              rot             = rotSiPM_mX;
              logicGreaseFace = logicGrease_XZ;
            } else if (face == 2) {//Y+
              greaseCenter    = G4ThreeVector(posX, posY + crystal_ly/2 + grease_t/2, posZ);
              faceCenter      = G4ThreeVector(posX, posY + crystal_ly/2 + grease_t + sipm_t/2, posZ);
              rot             = rotSiPM_pY;
              logicGreaseFace = logicGrease_Y;
            } else if (face == 3) {//Y-
              greaseCenter    = G4ThreeVector(posX, posY - crystal_ly/2 - grease_t/2, posZ);
              faceCenter      = G4ThreeVector(posX, posY - crystal_ly/2 - grease_t - sipm_t/2, posZ);
              rot             = rotSiPM_mY;
              logicGreaseFace = logicGrease_Y;
            } else if (face == 4) {//Z+
              greaseCenter    = G4ThreeVector(posX, posY, posZ + crystal_l/2 + grease_t/2);
              faceCenter      = G4ThreeVector(posX, posY, posZ + crystal_l/2 + grease_t + sipm_t/2);
              rot             = rotSiPM_pZ;
              logicGreaseFace = logicGrease_XZ;
            } else {//Z-  face == 5
              greaseCenter    = G4ThreeVector(posX, posY, posZ - crystal_l/2 - grease_t/2);
              faceCenter      = G4ThreeVector(posX, posY, posZ - crystal_l/2 - grease_t - sipm_t/2);
              rot             = rotSiPM_mZ;
              logicGreaseFace = logicGrease_XZ;
            }

            // One grease slab covering the entire crystal face
            char greaseNameBuf[64];
            std::snprintf(greaseNameBuf, sizeof(greaseNameBuf), "Grease_c%d_f%d", crystalId, face);
            new G4PVPlacement(rot, greaseCenter, logicGreaseFace, G4String(greaseNameBuf), logicEnv, false, crystalId * 6 + face, checkOverlaps);

            // 4×4 SiPM pixels on this face
            for (G4int j = 0; j < SiPm_np; ++j) {
              for (G4int k = 0; k < SiPm_np; ++k) {
                G4double u = (j - 1.5) * pitch;
                G4double v = (k - 1.5) * pitch;
                G4ThreeVector offset;
                if (face == 0 || face == 1) offset = G4ThreeVector(0, u, v);
                else if (face == 2 || face == 3) offset = G4ThreeVector(u, 0, v);
                else offset = G4ThreeVector(u, v, 0);
                G4ThreeVector posSiPM = faceCenter + offset;
                G4int copyNoSiPM = crystalId * 6 * 16 + face * 16 + j * 4 + k;
                char sipmNameBuf[64];
                std::snprintf(sipmNameBuf, sizeof(sipmNameBuf), "SiPM_%d", copyNoSiPM);
                new G4PVPlacement(rot, posSiPM, logicSiPM, G4String(sipmNameBuf), logicEnv, false, copyNoSiPM, checkOverlaps);
              }
            }
          }
        }
      }
    }

    fScoringVolume = logicCrystal;
    fCrystal_gap = Crystal_gap / mm;
    fCrystal_nx = Crystal_nx;
    fCrystal_ny = Crystal_ny;
    fCrystal_nz = Crystal_nz;
    fcrystal_l = crystal_l / mm;
    fcrystal_ly = crystal_ly / mm;
    fCrystal_x = Crystal_x;
    fCrystal_y = Crystal_y;
    fCrystal_z = Crystal_z;


      //Crystal Optical Surface
      G4OpticalSurface* crystalsurface = new G4OpticalSurface("CrystalSurface");
      crystalsurface->SetType(dielectric_dielectric);
      crystalsurface->SetModel(unified); 
      crystalsurface->SetFinish(polished);
      // crystalsurface->SetFinish(ground);
      crystalsurface->SetSigmaAlpha(Surface_Sigma);
      new G4LogicalSkinSurface("CrystalSurface",logicCrystal, crystalsurface);

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

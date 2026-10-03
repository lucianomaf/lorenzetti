
#include "DetectorConstruction_v1.h"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4ProductionCuts.hh"
#include "G4AutoDelete.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UniformMagField.hh"
#include "G4ProductionCuts.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4RegionStore.hh"
#include <string>
#include <sstream>

namespace{
template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}
}

G4ThreadLocal G4GlobalMagFieldMessenger* DetectorConstruction_v1::m_magFieldMessenger = nullptr;


/**
 * @class DetectorConstruction_v1
 * @brief Defines the full detector geometry for the simulation.
 * 
 * This class implements the G4VUserDetectorConstruction interface to build the 
 * Lorenzetti detector (calorimeters) using Geant4 geometry primitives.
 * It manages materials, logical volumes, and placements for both horizontal 
 * and vertical calorimeter plates.
 * 
 * Properties:
 * - UseMagneticField: Toggle global magnetic field (2 Tesla).
 * - UseSolenoidField: Uniform 2 Tesla axial field confined to the volume inside the ATLAS
 *   central solenoid (inner radius 1.23 m, axial length 5.8 m), with no field elsewhere.
 *   Cannot be combined with UseMagneticField.
 * - CutOnPhi: Restrict the detector to a phi wedge (used for debugging/visualization).
 */
DetectorConstruction_v1::DetectorConstruction_v1(std::string name)
 : 
  IMsgService(name), 
   G4VUserDetectorConstruction(),
   m_checkOverlaps(true)
{
  declareProperty( "UseMagneticField"           , m_useMagneticField=true     );
  declareProperty( "UseSolenoidField"           , m_useSolenoidField=false    );
  declareProperty( "CutOnPhi"                 , m_cutOnPhi=false            );
  declareProperty( "OutputLevel"                , m_outputLevel=0             ); 
}


DetectorConstruction_v1::~DetectorConstruction_v1()
{;}



//
// Add volume
//
void DetectorConstruction_v1::AddVolume(std::string region,
                                        int plates,
                                        std::string absorberMaterial,
                                        std::string gapMaterial,
                                        int nofLayers,
                                        double absoThickness,
                                        double gapThickness,
                                        double layerClearance,
                                        double rMin,
                                        double rMax,
                                        double zSize,
                                        double x,
                                        double y,
                                        double z,
                                        // production cuts
                                        double electronCut,
                                        double positronCut,
                                        double gammaCut,
                                        double photonCut
                                        )
{
  m_volumes.push_back(Volume{region,plates,absorberMaterial,gapMaterial,nofLayers,absoThickness,gapThickness,layerClearance,rMin,rMax,zSize,x,y,z,
                      electronCut,positronCut,gammaCut,photonCut});
}

void DetectorConstruction_v1::SetAbsorberSplit(std::string region, double zSplit, std::string absorberMaterial2,
                                               double absoThickness2)
{
  m_absorberSplit[region] = AbsorberSplit{zSplit, absorberMaterial2, absoThickness2};
}



G4VPhysicalVolume* DetectorConstruction_v1::Construct()
{
  // Define materials
  DefineMaterials();
  // Define volumes

  G4Material* defaultMaterial  = G4Material::GetMaterial("Vacuum");

  //
  // World
  //
  G4VSolid* worldS = new G4Tubs( "World",     // its name
                                 0,           // R min
                                 6*m,         // R max
                                 10*m,        // Z max
                                 0*deg,       // phi_min
                                 (!m_cutOnPhi)?(360*deg):(235*deg)      // phi_max
                                 );

  G4LogicalVolume* worldLV = new G4LogicalVolume(
                                 worldS,           // its solid
                                 defaultMaterial,  // its material
                                 "World");         // its name

  G4VPhysicalVolume* worldPV = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 m_checkOverlaps);  // checking overlaps

  //
  // Solenoid field volume
  //
  // Vacuum cylinder filling the inside of the ATLAS central solenoid, which carries the
  // solenoid field (see ConstructSDandField). Dimensions from JINST 3 (2008) S08003, table 2.1:
  // inner diameter 2.46 m, axial length 5.8 m. The coil and its cryostat are not modelled.
  if (m_useSolenoidField){
    if (m_useMagneticField){
      MSG_FATAL("UseSolenoidField and UseMagneticField cannot be enabled together. Abort!");
    }
    const double solenoidInnerRadius = 1230*mm;
    const double solenoidHalfLength  = 2900*mm;
    G4VSolid* solenoidS = new G4Tubs( "SolenoidField",
                                      0,
                                      solenoidInnerRadius,
                                      solenoidHalfLength,
                                      0*deg,
                                      (!m_cutOnPhi)?(360*deg):(235*deg) );
    m_solenoidLV = new G4LogicalVolume( solenoidS, defaultMaterial, "SolenoidField" );
    new G4PVPlacement( 0, G4ThreeVector(), m_solenoidLV, "SolenoidField", worldLV, false, 0, m_checkOverlaps );
    MSG_INFO( "Creating solenoid field volume: R < " << solenoidInnerRadius/mm << " mm, |z| < "
              << solenoidHalfLength/mm << " mm" );
  }


  for (auto &volume : m_volumes){

    MSG_INFO( "Creating Volume with name " << volume.name );

    if(volume.plates == 0){ //Plates::Horizontal){

      if (volume.layerClearance != 0){
        MSG_FATAL("LayerClearance is only implemented for vertical plates (volume " << volume.name << "). Abort!");
      }
      auto region = GetRegion(volume.name);


      // Create a region
      CreateHorizontalPlates( worldLV,
                    volume.name,
                    defaultMaterial, // default
                    G4Material::GetMaterial(volume.absorberMaterial), // absorber
                    G4Material::GetMaterial(volume.gapMaterial), // gap
                    volume.nofLayers, // ATLAS-like
                    volume.absoThickness, // abso
                    volume.gapThickness, // gap
                    volume.rMin, // start radio,
                    volume.zSize ,// z
                    G4ThreeVector(volume.x,volume.y,volume.z),
                    region );
      G4ProductionCuts* cuts=new G4ProductionCuts();
      cuts->SetProductionCut(volume.gammaCut,"gamma");
      cuts->SetProductionCut(volume.electronCut,"e-");
      cuts->SetProductionCut(volume.positronCut,"e+");
      cuts->SetProductionCut(volume.photonCut,"proton");
      // assign cuts to the region and return succesfully
      region->SetProductionCuts(cuts);
      
    }else if (volume.plates == 1){ //Plates::Vertical){

      auto region = GetRegion(volume.name);
      // Create a region
      CreateVerticalPlates( worldLV,
                    volume.name,
                    defaultMaterial, // default
                    G4Material::GetMaterial(volume.absorberMaterial), // absorber
                    G4Material::GetMaterial(volume.gapMaterial), // gap
                    volume.nofLayers, // ATLAS-like
                    volume.absoThickness, // abso
                    volume.gapThickness, // gap
                    volume.layerClearance, // empty space per layer
                    volume.rMin, // start radio,
                    volume.rMax,
                    volume.zSize ,// z
                    G4ThreeVector(volume.x,volume.y,volume.z),
                    region );


      G4ProductionCuts* cuts=new G4ProductionCuts();
      cuts->SetProductionCut(volume.gammaCut,"gamma");
      cuts->SetProductionCut(volume.electronCut,"e-");
      cuts->SetProductionCut(volume.positronCut,"e+");
      cuts->SetProductionCut(volume.photonCut,"proton");
      // assign cuts to the region and return succesfully
      region->SetProductionCuts(cuts);
    
    }else{
      MSG_FATAL("Volume type is not recognize. Abort!");
    }
  }// Loop over volumes

  return worldPV;
}


void DetectorConstruction_v1::DefineMaterials()
{
  // Lead material defined using NIST Manager
  G4NistManager* nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Cu");
  nistManager->FindOrBuildMaterial("G4_Pb");
  nistManager->FindOrBuildMaterial("G4_Fe");
  nistManager->FindOrBuildMaterial("G4_Al");
  nistManager->FindOrBuildMaterial("G4_Si");
  nistManager->FindOrBuildMaterial("G4_CESIUM_IODIDE");
  nistManager->FindOrBuildMaterial("PLASTIC SCINTILLATOR");
  G4Element* elH = new G4Element("Hydrogen", "H", 1., 1.0079 * g/mole);
  G4Element* elC = new G4Element("Carbon", "C", 6., 12.011 * g/mole);
  G4Material* scinti = new G4Material("PLASTIC SCINTILLATOR", 1.032 * g/cm3, 2); // Organic plastic scintillator, polystyrene
  scinti->AddElement(elH, 8.5*perCent);
  scinti->AddElement(elC, 91.5*perCent);

  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density;
  // The argon by NIST Manager is a gas with a different density
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);

  // Vacuum
  new G4Material("Vacuum", z=1., a=1.01*g/mole,density= universe_mean_density,kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Absorbers of the ATLAS barrel EM calorimeter (option AtlasEmb, geometry/python/v1/ECAL.py), as homogeneous mixtures
  // with the masses per unit area of one absorber: lead (1.53 mm below |eta| = 0.8, 1.13 mm above), two 0.2 mm sheets of
  // stainless steel (taken as iron) and 0.832 mm of polyimide for the glue and the readout electrode (JINST 3 (2008)
  // S08003, sec. 5.2.1; the polyimide thickness is not given there and is fitted to the X0 of fig. 5.1 at eta = 0).
  // Not used by the default geometry.
  {
    G4Material *pb = nistManager->FindOrBuildMaterial("G4_Pb");
    G4Material *fe = nistManager->FindOrBuildMaterial("G4_Fe");
    G4Material *kapton = nistManager->FindOrBuildMaterial("G4_KAPTON");
    const double tfe = 0.4, tkapton = 0.832; // mm
    for( double tpb : {1.53, 1.13} ){
      const double mpb = tpb*pb->GetDensity(), mfe = tfe*fe->GetDensity(), mk = tkapton*kapton->GetDensity();
      const double m = mpb + mfe + mk;
      const std::string matName = (tpb > 1.3) ? "ATLAS_EMB_ABSORBER_153" : "ATLAS_EMB_ABSORBER_113";
      G4Material *mix = new G4Material(matName, m/(tpb + tfe + tkapton), 3);
      mix->AddMaterial(pb, mpb/m);
      mix->AddMaterial(fe, mfe/m);
      mix->AddMaterial(kapton, mk/m);
    }
  }

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}







/**
 * @brief Constructs the horizontal plates of a calorimeter module (Barrel).
 * 
 * Creates a "sandwich" structure of absorber, gap, and absorber layers arranged 
 * radially (for cylindrical/barrel geometry).
 * 
 * @param worldLV Mother logical volume.
 * @param name Unique name for the volume.
 * @param defaultMaterial Envelope material (usually Vacuum).
 * @param absorberMaterial Absorber material (e.g., Pb, Fe).
 * @param gapMaterial Active material (e.g., LAr, Scintillator).
 * @param nofLayers Number of repeating layers.
 * @param absoThickness Thickness of the absorber plate.
 * @param gapThickness Thickness of the active gap.
 * @param calorRmin Inner radius.
 * @param calorZ Total length in Z.
 * @param center_pos Center position (displacement).
 * @param region G4Region associated with this detector part.
 */
void DetectorConstruction_v1::CreateHorizontalPlates(  G4LogicalVolume *worldLV, 
                                                    std::string name,  
                                                    G4Material *defaultMaterial,
                                                    G4Material *absorberMaterial,
                                                    G4Material *gapMaterial,
                                                    int nofLayers,
                                                    double absoThickness,
                                                    double gapThickness,
                                                    double calorRmin,
                                                    double calorZ,
                                                    const G4ThreeVector &center_pos,
                                                    G4Region *region
                                                    ) 

{
  if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("DetectorConstruction_v1::DefineVolumes()", "MyCode0001", FatalException, msg);
  }


  G4double layerThickness=absoThickness+gapThickness; 

  G4VSolid* calorimeterS = new G4Tubs( name,// its name
                                 calorRmin, // R min 1700mm
                                 calorRmin+ nofLayers*layerThickness, // R max 48cm+1700mm
                                 calorZ/2,    // Z max, +250cm
                                 0*deg,     // phi_min
                                 (!m_cutOnPhi)?(360*deg):(235*deg)    // phi_max
                                 ) ;




  G4LogicalVolume* calorLV = new G4LogicalVolume( calorimeterS,     // its solid
                                                  defaultMaterial,  // its material
                                                  name);   // its name
  new G4PVPlacement(
                 0,                 // no rotation
                 center_pos,        // at (0,0,0)
                 calorLV,           // its logical volume
                 name,              // its name
                 worldLV,           // its mother  volume
                 false,             // no boolean operation
                 0,                 // copy number
                 m_checkOverlaps);  // checking overlaps
  
  region->AddRootLogicalVolume(calorLV);



  auto split = m_absorberSplit.find(name);
  for (G4int layer=0; layer < nofLayers; ++layer){

    G4VSolid* layerS = new G4Tubs(name+ "_Layer",// its name
                                 calorRmin + layer*layerThickness,        // R min 1700mm
                                 calorRmin + (layer+1)*layerThickness,   // R max 48cm+1700mm
                                 calorZ/2,            // Z max, +250cm
                                 0*deg,             // phi_min
                                 (!m_cutOnPhi)?(360*deg):(235*deg)            // phi_max
                                 ) ;

    G4LogicalVolume* layerLV = new G4LogicalVolume(
                                                  layerS,           // its solid
                                                  defaultMaterial,  // its material
                                                  name+"_Layer");   // its name


    new G4PVPlacement(
                 0,                 // no rotation
                 G4ThreeVector(0,0,0),
                 //center_pos,        // at (0,0,0)
                 layerLV,           // its logical volume
                 name+"_Layer",     // its name
                 calorLV,           // its mother  volume
                 false,             // no boolean operation
                 0,                 // copy number
                 m_checkOverlaps);  // checking overlaps
    


    if( split == m_absorberSplit.end() ){
    G4VSolid* absorverS = new G4Tubs( name+"_Abso",// its name
                                 calorRmin + layer*(absoThickness + gapThickness) ,     // R min 1700mm
                                 calorRmin + layer*(absoThickness + gapThickness) + absoThickness,   // R max 48cm+1700mm
                                 calorZ/2,          // Z max, +250cm
                                 0*deg,                 // phi_min
                                 (!m_cutOnPhi)?(360*deg):(235*deg)             // phi_max
                                 ) ;

    G4LogicalVolume* absorverLV = new G4LogicalVolume(
                                                  absorverS,           // its solid
                                                  absorberMaterial,  // its material
                                                  name+"_Abso");         // its name


    new G4PVPlacement(
                 0,                 // no rotation
                 G4ThreeVector(0,0,0),
                 //center_pos,        // at (0,0,0)
                 absorverLV,        // its logical volume
                 name+"_Abso",      // its name
                 layerLV,           // its mother  volume
                 false,             // no boolean operation
                 0,                 // copy number
                 m_checkOverlaps);  // checking overlaps




    G4VSolid* gapS = new G4Tubs( name+"_Gap",// its name
                                 calorRmin + layer*(absoThickness + gapThickness) + absoThickness,   // R max 48cm+1700mm
                                 calorRmin + (layer+1)*(absoThickness + gapThickness),   // R max 48cm+1700mm
                                 calorZ/2,          // Z max, +250cm
                                 0*deg,           // phi_min
                                 (!m_cutOnPhi)?(360*deg):(235*deg)          // phi_max
        ); 

    G4LogicalVolume* gapLV = new G4LogicalVolume(
                                                  gapS,           // its solid
                                                  gapMaterial,  // its material
                                                  name+"_Gap");         // its name


    new G4PVPlacement(
                 0,                 // no rotation
                 G4ThreeVector(0,0,0),
                 //center_pos,        // at (0,0,0)
                 gapLV,             // its logical volume
                 name+"_Gap",       // its name
                 layerLV,           // its mother  volume
                 false,             // no boolean operation
                 0,                 // copy number
                 m_checkOverlaps);  // checking overlaps


    } else {
      // Absorber changing for |z| > zSplit (SetAbsorberSplit): three pieces in z, the same layer thickness, the gap
      // taking the rest of the layer in each piece.
      const AbsorberSplit &sp = split->second;
      G4Material *absorberMaterial2 = G4Material::GetMaterial(sp.material);
      const double zh = calorZ/2;
      if( !absorberMaterial2 || sp.z <= 0 || sp.z >= zh || sp.thickness <= 0 || sp.thickness >= layerThickness ){
        G4ExceptionDescription msg;
        msg << "Invalid absorber split for " << name;
        G4Exception("DetectorConstruction_v1::CreateHorizontalPlates()", "MyCode0002", FatalException, msg);
      }
      const double r0 = calorRmin + layer*layerThickness;
      struct Piece { double halfz; double zc; G4Material *mat; double t; };
      const Piece pieces[3] = { {sp.z, 0., absorberMaterial, absoThickness},
                                {(zh - sp.z)/2,  (zh + sp.z)/2, absorberMaterial2, sp.thickness},
                                {(zh - sp.z)/2, -(zh + sp.z)/2, absorberMaterial2, sp.thickness} };
      for( const auto &piece : pieces ){
        G4VSolid* aS = new G4Tubs( name+"_Abso", r0, r0 + piece.t, piece.halfz, 0*deg, (!m_cutOnPhi)?(360*deg):(235*deg) );
        G4LogicalVolume* aLV = new G4LogicalVolume( aS, piece.mat, name+"_Abso" );
        new G4PVPlacement( 0, G4ThreeVector(0,0,piece.zc), aLV, name+"_Abso", layerLV, false, 0, m_checkOverlaps );
        G4VSolid* gS = new G4Tubs( name+"_Gap", r0 + piece.t, r0 + layerThickness, piece.halfz, 0*deg, (!m_cutOnPhi)?(360*deg):(235*deg) );
        G4LogicalVolume* gLV = new G4LogicalVolume( gS, gapMaterial, name+"_Gap" );
        new G4PVPlacement( 0, G4ThreeVector(0,0,piece.zc), gLV, name+"_Gap", layerLV, false, 0, m_checkOverlaps );
      }
    }
  }// Loop over calorimeter layers
}


void DetectorConstruction_v1::CreateVerticalPlates(  G4LogicalVolume *worldLV, 
                                                  std::string name,  
                                                  G4Material *defaultMaterial,
                                                  G4Material *absorberMaterial,
                                                  G4Material *gapMaterial,
                                                  int nofLayers,
                                                  double absoThickness,
                                                  double gapThickness,
                                                  double layerClearance,
                                                  double calorRmin,
                                                  double calorRmax,
                                                  double calorZ,
                                                  const G4ThreeVector &center_pos,
                                                  G4Region *region
                                                  ) 
{
  if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("DetectorConstruction_v1::DefineVolumes()", "MyCode0001", FatalException, msg);
  }

  const bool evenNofLayers = !(nofLayers % 2);
  const bool invertOrder = center_pos.z() < 0.;

  // Each layer is absorber + gap + clearance. The absorber sits at one end of the layer and the
  // gap at the other, so a non-zero clearance leaves an empty space (envelope material) between
  // the gap and the absorber of the same layer: along the stacking axis the sequence is
  // absorber, gap, clearance, absorber, ... With clearance = 0 the layout is unchanged.
  // Used by the ATLAS-like tile calorimeter: 14 mm steel + 3 mm tile + 1 mm left by the 3 mm
  // tile in its 4 mm pocket, period 18 mm (TILECAL, hep-ex/9904032; JINST 3 S08003, sec. 5.3.1.2).
  G4double layerThickness=absoThickness+gapThickness+layerClearance; 

  G4VSolid* calorimeterS = new G4Tubs( name,// its name
                                 calorRmin,
                                 calorRmax,
                                 calorZ/2, // Size in z
                                 0*deg,     // phi_min
                                 (!m_cutOnPhi)?(360*deg):(235*deg)    // phi_max
                                 ) ;

  G4LogicalVolume* calorLV = new G4LogicalVolume( calorimeterS,     // its solid
                                                  defaultMaterial,  // its material
                                                  name);   // its name


  new G4PVPlacement(
                 0,                 // no rotation
                 center_pos,        // at (0,0,0)
                 calorLV,           // its logical volume
                 name,              // its name
                 worldLV,           // its mother  volume
                 false,             // no boolean operation
                 0,                 // copy number
                 m_checkOverlaps);  // checking overlaps
  
  region->AddRootLogicalVolume(calorLV);

  for (G4int layer=-nofLayers/2; layer <= nofLayers/2; ++layer){
    if ( layer == 0 && evenNofLayers){
      continue;
    }

    G4VSolid* layerS = new G4Tubs(name+ "_Layer",// its name
                                 calorRmin,
                                 calorRmax,
                                 layerThickness/2,
                                 0*deg,
                                 (!m_cutOnPhi)?(360*deg):(235*deg)
                                 );

    G4LogicalVolume* layerLV = new G4LogicalVolume(
                                                  layerS,           // its solid
                                                  defaultMaterial,  // its material
                                                  name+"_Layer");   // its name

    const auto layerLVCenter = G4ThreeVector(0,0,
                   ((invertOrder)?-1.:+1.) * 
                   (sign(-layer)*((evenNofLayers)?0.5:0.) + (layer))
                   *layerThickness // center at with respect to the calorimeter center
                 );

    //G4cout << "Creating layer " << layer << " for: " << layerLVCenter << ". Size: " << layerThickness << std::endl;

    new G4PVPlacement(
                 0,                 // no rotation
                 layerLVCenter,     // center
                 layerLV,           // its logical volume
                 name+"_Layer",     // its name
                 calorLV,           // its mother  volume
                 false,             // no boolean operation
                 0,                 // copy number
                 m_checkOverlaps);  // checking overlaps



    G4VSolid* absorverS = new G4Tubs( name+"_Abso",// its name
                                 calorRmin,     // R min
                                 calorRmax,     // R max
                                 absoThickness/2,  // 
                                 0*deg,         // phi_min
                                 (!m_cutOnPhi)?(360*deg):(235*deg)        // phi_max
                                 ) ;

    G4LogicalVolume* absorverLV = new G4LogicalVolume(
                                                  absorverS,           // its solid
                                                  absorberMaterial,  // its material
                                                  name+"_Abso");         // its name


    new G4PVPlacement(
                 0,                 // no rotation
                 G4ThreeVector(0,0,
                   ((invertOrder)?-1.:+1) * 
                   layerThickness/2
                   + ((invertOrder)?+1.:-1)/2 * 
                    absoThickness
                 ),
                 absorverLV,        // its logical volume
                 name+"_Abso",      // its name
                 layerLV,           // its mother  volume
                 false,             // no boolean operation
                 0,                 // copy number
                 m_checkOverlaps);  // checking overlaps




    G4VSolid* gapS = new G4Tubs( name+"_Gap",// its name
                                 calorRmin,   // R max
                                 calorRmax,   // R max
                                 gapThickness/2,   // Z max,
                                 0*deg,           // phi_min
                                 (!m_cutOnPhi)?(360*deg):(235*deg)          // phi_max
        ); 

    G4LogicalVolume* gapLV = new G4LogicalVolume(
                                                  gapS,           // its solid
                                                  gapMaterial,  // its material
                                                  name+"_Gap");         // its name


    new G4PVPlacement(
                 0,                 // no rotation
                 G4ThreeVector(0,0,
                   ((invertOrder)?+1.:-1) * 
                   layerThickness/2
                   + ((invertOrder)?-1.:+1) * 
                    gapThickness/2
                 ),
                 gapLV,             // its logical volume
                 name+"_Gap",       // its name
                 layerLV,           // its mother  volume
                 false,             // no boolean operation
                 0,                 // copy number
                 m_checkOverlaps);  // checking overlaps
  }// Loop over calorimeter layers

}





void DetectorConstruction_v1::ConstructSDandField(){

  if (m_useMagneticField){
    MSG_INFO("Set magnetic field")
    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    G4ThreeVector fieldValue(0.0 , 0.0, 2*tesla);
    m_magFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
    m_magFieldMessenger->SetVerboseLevel(1);

    // Register the field messenger for deleting
    G4AutoDelete::Register(m_magFieldMessenger);
  }

  if (m_useSolenoidField && m_solenoidLV){
    MSG_INFO("Set solenoid magnetic field (2 T, only inside the solenoid volume)")
    // Local field: attached to the solenoid volume only, so particles leaving it
    // (and the showers in the calorimeters) are not bent.
    auto solenoidField = new G4UniformMagField( G4ThreeVector(0.0, 0.0, 2*tesla) );
    auto fieldManager  = new G4FieldManager();
    fieldManager->SetDetectorField( solenoidField );
    fieldManager->CreateChordFinder( solenoidField );
    m_solenoidLV->SetFieldManager( fieldManager, true );
    // The field manager is owned and deleted by the G4FieldManagerStore: registering it
    // with G4AutoDelete as well deletes it twice at exit.
    G4AutoDelete::Register( solenoidField );
  }

}



G4Region* DetectorConstruction_v1::GetRegion( std::string name ){
  auto reg = G4RegionStore::GetInstance()->GetRegion(name);
  return reg ? reg : new G4Region(name);
}
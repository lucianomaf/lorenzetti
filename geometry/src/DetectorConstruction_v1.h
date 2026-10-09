#ifndef DetectorConstruction_v1_h
#define DetectorConstruction_v1_h

#include "GaugiKernel/MsgStream.h"
#include "GaugiKernel/Property.h"
#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4Region.hh"
#include "globals.hh"
#include "G4Cache.hh"
#include <map>
#include <set>
#include <vector>


//enum Plates{
//    Horizontal       = 0,
//    Vertical         = 1,
//};

class G4VPhysicalVolume;

class G4GlobalMagFieldMessenger;

class DetectorConstruction_v1 : public G4VUserDetectorConstruction, public MsgService, public Gaugi::PropertyService
{
  public:


    struct Volume{
      std::string name;
      int plates;
      std::string absorberMaterial;
      std::string gapMaterial;
      int nofLayers;
      double absoThickness;
      double gapThickness;
      double layerClearance; // empty space (envelope material) per layer, vertical plates only
      double rMin;
      double rMax;
      double zSize;
      double x;
      double y;
      double z;

      double electronCut;
      double positronCut;
      double gammaCut;
      double photonCut;
    };

    DetectorConstruction_v1(std::string);
    virtual ~DetectorConstruction_v1();
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    // get methods
    const G4VPhysicalVolume* GetAbsorberPV() const;
    const G4VPhysicalVolume* GetGapPV() const;

    void AddVolume(std::string region,
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
                   );

    // Horizontal plates whose absorber changes for |z| > zSplit, with the same layer thickness (the gap takes the rest
    // of the layer): used for the lead of the ATLAS barrel EM calorimeter (1.53 mm below |eta| = 0.8, 1.13 mm above;
    // option AtlasEmb in geometry/python). Must be called after AddVolume for the same region.
    void SetAbsorberSplit(std::string region, double zSplit, std::string absorberMaterial2, double absoThickness2);
    struct AbsorberSplit{ double z; std::string material; double thickness; };

    // A material from elements, defined by the geometry in Python (density in g/cm3, mass fractions by element symbol):
    // the phi-averaged mixtures of the ATLAS pieces in front of the barrel accordion (F41 conserto 17c).
    void AddMixture(std::string name, double density, std::vector<std::string> symbols, std::vector<double> fractions);
    // Plates::Polycone (2): one G4Polycone of the absorber material with these planes. Call after AddVolume.
    void SetPolycone(std::string region, std::vector<double> z, std::vector<double> rmin, std::vector<double> rmax);
    // Place this volume inside the solenoid field volume (UseSolenoidField) instead of the world. Call after AddVolume.
    void SetInSolenoidField(std::string region);
    struct Mixture{ std::string name; double density; std::vector<std::string> symbols; std::vector<double> fractions; };
    struct PolyconePlanes{ std::vector<double> z, rmin, rmax; };

  private:

    std::vector<Volume> m_volumes;
    std::map<std::string, AbsorberSplit> m_absorberSplit;
    std::vector<Mixture> m_mixtures;
    std::map<std::string, PolyconePlanes> m_polycones;
    std::set<std::string> m_inSolenoidField;
    double m_solenoidFieldRadius;

    void CreatePolycone( G4LogicalVolume *motherLV, std::string name, G4Material *material, const PolyconePlanes &planes,
                         G4Region *region );

    // methods
    void DefineMaterials();
    
    G4VPhysicalVolume* DefineVolumes();


    void CreateHorizontalPlates(  G4LogicalVolume *worldLV, 
                                  std::string name,  
                                  G4Material *defaltMaterial,
                                  G4Material *absorberMaterial,
                                  G4Material *gapMaterial,
                                  int nofLayers,
                                  double absoThickness,
                                  double gapThickness,
                                  double calorRmin,
                                  double calorZ,
                                  const G4ThreeVector &center_pos,
                                  G4Region* region);

    void CreateVerticalPlates(  G4LogicalVolume *worldLV, 
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
                                G4Region *region);

    G4Region* GetRegion( std::string name );


    bool m_checkOverlaps; // option to activate checking of volumes overlaps
    bool m_useMagneticField;
    bool m_useSolenoidField;
    bool m_cutOnPhi;
    int m_outputLevel;

    // Field volume of the solenoid (only built with UseSolenoidField)
    G4LogicalVolume* m_solenoidLV = nullptr;

    static G4ThreadLocal G4GlobalMagFieldMessenger*  m_magFieldMessenger;
};




#endif


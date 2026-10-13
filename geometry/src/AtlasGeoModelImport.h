#ifndef AtlasGeoModelImport_h
#define AtlasGeoModelImport_h

/**
 * Volumes of the ATLAS geometry of the simulation read from its GeoModel SQLite dump and converted to Geant4 (GeoModelRead and
 * GeoModel2G4 of the GeoModel project), for --atlas-run4 (F41 conserto 27): the inner detector of Run 4 (ITk, HGTD), the beam pipe and
 * the end-cap cryostats of the LAr calorimeter.
 *
 * Built only with the GeoModel libraries (CMake: GEOMODEL_DIR; LORENZETTI_WITH_GEOMODEL). The header has no GeoModel type, so that the
 * rest of the code is compiled without them.
 */

#include "G4Transform3D.hh"
#include <memory>
#include <set>
#include <string>
#include <vector>

class G4LogicalVolume;

class AtlasGeoModelImport
{
  public:
    // A cloned ATLAS volume (one logical volume per placement) that can hold Lorenzetti volumes: its global transformation and depth
    struct Container {
      G4LogicalVolume *lv;
      G4Transform3D x;
      int depth;
      std::string name;
    };

    explicit AtlasGeoModelImport( const std::string &dbFile );
    ~AtlasGeoModelImport();

    // true when Lorenzetti was built with the GeoModel libraries
    static bool available();

    // The children of the world with these logical volume names, whole, placed in mother with their transformation.
    void placeWorldChildren( G4LogicalVolume *mother, const std::set<std::string> &names, bool checkOverlaps );

    // The children of the world volumes named `name` placed in mother without the volume itself (its own material is the air of the
    // world), with the composed transformation. The descendants whose logical volume name starts with one of `dropped` are left out
    // (Lorenzetti builds them); every volume that has one of them below it is cloned per placement and returned as a container.
    std::vector<Container> placeUnwrapped( G4LogicalVolume *mother, const std::string &name, const std::vector<std::string> &dropped,
                                           bool checkOverlaps );

    // Logical volumes converted so far whose name (with the prefix of their subdetector) is in `names`: for the regions of ATLAS
    std::vector<G4LogicalVolume*> logicalVolumes( const std::set<std::string> &names ) const;

  private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

#endif

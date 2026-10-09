#ifndef CaloHitMaker_h
#define CaloHitMaker_h

#include "GaugiKernel/Algorithm.h"
#include "GaugiKernel/DataHandle.h"
#include "GaugiKernel/AlgTool.h"
#include "CaloHit/CaloHit.h"



/**
 * @class CaloHitMaker
 * @brief Algorithm to build calorimeter hits during simulation.
 * 
 * Collects energy deposits from Geant4 steps and integrates them into hits
 * corresponding to readout cells.
 */
class CaloHitMaker : public Gaugi::Algorithm
{

  typedef std::map<unsigned long int, xAOD::CaloHit* > collection_map_t;

  public:
  
    /** Contructor **/
    CaloHitMaker( std::string name );
    /** Destructor **/
    ~CaloHitMaker()=default;
    
    /*! initialize the algorithm **/
    virtual StatusCode initialize() override;
    
    /*! Book all histograms into the current storegate **/
    virtual StatusCode bookHistograms( SG::EventContext &ctx ) const override;
    
    /*! Execute in step action step from geant core **/
    virtual StatusCode execute( SG::EventContext &ctx , const G4Step *step) const override;
    
    /*! Execute in ComponentAccumulator **/
    virtual StatusCode execute( SG::EventContext &ctx , int /*evt*/ ) const override;
    
    /*! execute before start the step action **/
    virtual StatusCode pre_execute( SG::EventContext &ctx ) const override;
    
    /*! execute after the step action **/ 
    virtual StatusCode post_execute( SG::EventContext &ctx ) const override;
    
    /*! fill histogram in the end **/
    virtual StatusCode fillHistograms( SG::EventContext &ctx ) const override;
    
    virtual StatusCode finalize() override;

    /*! ATLAS-like EM end-cap (EmecCompartment > 0; geometry/python/v1/EMEC.py, --atlas-emec): the compartment (1-11)
        of the point (x, y, z), in mm, as in the ATLAS simulation (LArG4EC EnergyCalculator::FindIdentifier_Default),
        and, when it is the compartment of this maker, the eta index and the phi bin of its cell; 0 when the point is
        outside the wheel of this maker. Public so that the assignment can be tested from Python. */
    int emecFindCell( double x, double y, double z, int &etaBin, int &phiBin ) const;

    /*! ATLAS-like barrel EM and presampler cells (EmbMode > 0; geometry/python/v1/ECAL.py, --atlas-emb-cells): the cell
        (index in the table CellEta) of the point (x, y, z), in mm, for the sampling and region of this maker, or -1 when
        the point has no cell of that region (see the .cxx). Public so that the assignment can be tested from Python. */
    int embFindCell( double x, double y, double z ) const;

  private:
   
    int find( const std::vector<float> &vec, float value) const;
    unsigned long int hash(unsigned bin) const;

    /*! collection key */
    std::string m_collectionKey; // output

    std::vector<float> m_etaBins; 
    std::vector<float> m_phiBins; 
    float m_rMin;
    float m_rMax; 
    float m_zMin;
    float m_zMax;
    float m_noiseStd;
    // Keep only the energy deposited in the active medium (liquid argon, scintillator)
    bool m_activeEnergyOnly;
    // Apply Birks' law to the energy of each step (scintillator and liquid argon), as in ATLAS
    bool m_birksLaw;
    float birksEnergy( const G4Step *step ) const;

    /*! Sampling id for this reconstruction */
    int m_sampling;
    /*! Sampling id segment */
    int m_segment;
    /*! Detector id for this sampling*/
    int m_detector;
    /*! The start bunch crossing id for energy estimation */
    int m_bcid_start;
    /*! The end bunch crossing id for energy estimation */
    int m_bcid_end;
    /*! The time space (in ns) between two bunch crossings */
    float m_bc_duration;
    /*! Base histogram path */
    std::string m_histPath;
    /*! detailed histogram flags */
    bool m_detailedHistograms;

    unsigned int m_nEtaBins;
    unsigned int m_nPhiBins;

    // Cells given as (r, z) boxes instead of the eta x phi grid (ATLAS tile cells; empty by default).
    // One entry per cell: eta of its centre and its delta eta. One entry per box: r and z limits and the
    // index of the cell it belongs to (a cell may be made of several boxes, as BC in the tile barrel).
    std::vector<float> m_cellEta;
    std::vector<float> m_cellDeltaEta;
    std::vector<float> m_cellBoxRMin;
    std::vector<float> m_cellBoxRMax;
    std::vector<float> m_cellBoxZMin;
    std::vector<float> m_cellBoxZMax;
    std::vector<int>   m_cellBoxIndex;
    bool useCells() const { return !m_cellEta.empty(); }

    // Dual readout of the tile cells, as in ATLAS (TileGeoG4SDCalc::MakePmtEdepTime): each cell is read by two
    // PMTs, one on each side of the tiles in phi, and the energy of each step is shared between them according
    // to its azimuthal position across the module. 0 = off (default), 1 = ATLAS U-shape, 2 = linear sharing.
    int m_tileDualReadout;
    bool isTile() const { return m_sampling >= 5 && m_sampling <= 10; }
    void tilePmtWeights( float phiLocal, float z, float &w0, float &w1 ) const;
    int findCell( float radius, float z ) const;
    // Radius used to find the (r, z) box of a step: r itself (CellRadiusModules = 0, the default) or, with N modules,
    // the radius on the axis of the module that contains phi, r cos(phi - phi_c), modules of 2 pi / N from phi = 0
    // (moduleY of the ATLAS HEC simulation, LArG4HEC HECGeometry::CalculateIdentifier; ATLAS-like HEC, --atlas-hec).
    int m_cellRadiusModules;
    float cellRadius( float radius, float phi ) const;

    // ATLAS-like EM end-cap (see emecFindCell): the compartment read by this maker (1-11 as in the table s_geometry of
    // EnergyCalculator.cc; 0 = off, the default), its eta scale, offset and last index, the distance from the
    // mechanical to the electric focal point, the boundaries in z between the compartments (EmecSamplingSep: ZSEP12,
    // ZSEP23, ZIW, from the electric focal point, integers in 1e-4 cm as in the database) and the volumes (r, z boxes)
    // of the wheel of the compartment.
    int m_emecCompartment;
    float m_emecEtaScale;
    float m_emecEtaOffset;
    int m_emecMaxEta;
    float m_emecFocalShift;
    std::vector<int> m_emecZSep12;
    std::vector<int> m_emecZSep23;
    std::vector<int> m_emecZInner;
    // z boundary (mm) from a table value in 1e-4 cm, as the ATLAS simulation reads it (value in cm times CLHEP::cm)
    static double emecZ( int value ) { return ( value / 1e4 ) * 10.; }
    std::vector<float> m_emecWheelRMin;
    std::vector<float> m_emecWheelRMax;
    std::vector<float> m_emecWheelZMin;
    std::vector<float> m_emecWheelZMax;

    // ATLAS-like barrel EM and presampler cells (see embFindCell): mode (0 = off, the default; 1 = barrel; 2 = end-cap
    // presampler), sampling (0-3) and region (0 or 1) of the ATLAS identifier dictionary read by this maker.
    int m_embMode;
    int m_embSampling;
    int m_embRegion;

    // Free-running hits (simu_trf.py --free-running-hits): every deposit of the ATLAS tile cells also goes, per PMT, with
    // the identifier and the time of the ATLAS simulation, to a per-event list written in the layout of the ATLAS HITS
    // ntuple (CaloFreeRunningHitWriter). See freeRunningTile.
    bool m_freeRunningHits;
    std::string m_freeRunningListKey;
    bool m_freeRunningBinning;
    // ATLAS identifier fields and (r, z) centre of each cell (same order as m_cellEta; geometry/python/v1/TILE.py)
    std::vector<int>   m_cellAtlasSection;
    std::vector<int>   m_cellAtlasTower;
    std::vector<int>   m_cellAtlasSampling;
    std::vector<float> m_cellRCentre;
    std::vector<float> m_cellZCentre;
    void freeRunningTile( SG::EventContext &ctx, const G4Step *step, int cell, int phiBin, float edep, float w0,
                          float w1 ) const;
};


#endif

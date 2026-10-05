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



#include "CaloHit/CaloHit.h"
#include "CaloHit/CaloHitCollection.h"
#include "EventInfo/EventInfoContainer.h"
#include "CaloHitMaker.h"
#include "CaloFreeRunningHitList.h"
#include "TileUShape.h"
#include <algorithm>

#include "G4Kernel/CaloPhiRange.h"
#include "G4Kernel/constants.h"
#include "TVector3.h"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4Material.hh"
#include <cmath>

#include "TH1F.h"
#include "TH2F.h"
#include "TH2Poly.h"
#include "TGraph.h"
#include <cstdlib>

using namespace Gaugi;
using namespace SG;




/**
 * @class CaloHitMaker
 * @brief Collects energy deposits (Hits) during Geant4 simulation.
 * 
 * This algorithm runs at every Geant4 step. It checks if the energy deposit
 * occurred within the defined volume (Sampling) and integrates the energy
 * into a CaloHit object. It maps geometric position (x,y,z) to readout 
 * identifiers (eta, phi bins).
 * 
 * Properties:
 * - OutputCollectionKey: StoreGate key for the collection of hits.
 * - Eta/PhiBins: Readout segmentation.
 * - RMin/Max, ZMin/Max: Spatial boundaries of the sensitive volume.
 * - Sampling: Identifier for the calorimeter layer.
 */
CaloHitMaker::CaloHitMaker( std::string name ) : 
  IMsgService(name),
  Algorithm()
{
  declareProperty( "OutputCollectionKey"      , m_collectionKey="CaloHitCollection"   ); // output
  
  declareProperty( "EtaBins"                  , m_etaBins                             );
  declareProperty( "PhiBins"                  , m_phiBins                             );
  declareProperty( "RMin"                     , m_rMin                                );
  declareProperty( "RMax"                     , m_rMax                                );
  declareProperty( "ZMin"                     , m_zMin                                );
  declareProperty( "ZMax"                     , m_zMax                                );
  declareProperty( "Sampling"                 , m_sampling                            );
  declareProperty( "Segment"                  , m_segment                             );
  declareProperty( "Detector"                 , m_detector                            );
  declareProperty( "BunchIdStart"             , m_bcid_start=-7                       );
  declareProperty( "BunchIdEnd"               , m_bcid_end=8                          );
  declareProperty( "BunchDuration"            , m_bc_duration=25                      );
  declareProperty( "OutputLevel"              , m_outputLevel=1                       );
  declareProperty( "DetailedHistograms"       , m_detailedHistograms=false            );
  declareProperty( "HistogramPath"            , m_histPath="/CaloHitMaker"            );
  declareProperty( "SamplingNoiseStd"         , m_noiseStd=0                          );
  // When true, only steps in the active medium (liquid argon or plastic scintillator) are kept,
  // as in the ATLAS hits; by default the whole energy deposited in the cell volume is kept.
  declareProperty( "ActiveEnergyOnly"         , m_activeEnergyOnly=false              );
  // When true, Birks' law is applied to the energy of each step in the scintillator and in the
  // liquid argon, with the ATLAS simulation constants (see birksEnergy).
  declareProperty( "BirksLaw"                 , m_birksLaw=false                      );
  // Cells as (r, z) boxes (ATLAS tile cells, set from geometry/python/v1/TILE.py). When CellEta is empty
  // (the default) the cells are the eta x phi grid given by EtaBins and PhiBins.
  declareProperty( "CellEta"                  , m_cellEta                             );
  declareProperty( "CellDeltaEta"             , m_cellDeltaEta                        );
  declareProperty( "CellBoxRMin"              , m_cellBoxRMin                         );
  declareProperty( "CellBoxRMax"              , m_cellBoxRMax                         );
  declareProperty( "CellBoxZMin"              , m_cellBoxZMin                         );
  declareProperty( "CellBoxZMax"              , m_cellBoxZMax                         );
  declareProperty( "CellBoxIndex"             , m_cellBoxIndex                        );
  // With N > 0, the box of a step is found with the radius on the axis of its module (N modules from phi = 0),
  // r cos(phi - phi_c), as the ATLAS HEC simulation does (see cellRadius); 0 (the default) uses r.
  declareProperty( "CellRadiusModules"        , m_cellRadiusModules=0                 );
  // Dual readout of the tile cells (two PMTs per cell, as in ATLAS): 0 = off, 1 = ATLAS U-shape, 2 = linear.
  // Only for the tile samplings and only with the ATLAS tile cells (CellEta not empty); see tilePmtWeights.
  declareProperty( "TileDualReadout"          , m_tileDualReadout=0                   );
  // Free-running hits (simu_trf.py --free-running-hits): the deposits of the ATLAS tile cells also go, per PMT and with
  // the identifier and the time of the ATLAS simulation, to the list FreeRunningListKey (see freeRunningTile).
  declareProperty( "FreeRunningHits"          , m_freeRunningHits=false               );
  declareProperty( "FreeRunningListKey"       , m_freeRunningListKey="FreeRunningHits");
  declareProperty( "FreeRunningBinning"       , m_freeRunningBinning=true             );
  declareProperty( "CellAtlasSection"         , m_cellAtlasSection                    );
  declareProperty( "CellAtlasTower"           , m_cellAtlasTower                      );
  declareProperty( "CellAtlasSampling"        , m_cellAtlasSampling                   );
  declareProperty( "CellRCentre"              , m_cellRCentre                         );
  declareProperty( "CellZCentre"              , m_cellZCentre                         );



}

//!=====================================================================

StatusCode CaloHitMaker::initialize()
{
  MSG_INFO("1 Initialize..." << m_outputLevel);

  CHECK_INIT();
  setMsgLevel( (MSG::Level)m_outputLevel );
  m_nEtaBins = m_etaBins.size() - 1;
  m_nPhiBins = m_phiBins.size() - 1;
  if( m_tileDualReadout > 0 && isTile() && !useCells() ){
    MSG_FATAL( "TileDualReadout needs the ATLAS tile cells (CellEta)." );
  }
  if( m_freeRunningHits && isTile() && ( !useCells() || m_tileDualReadout <= 0 ||
                                        m_cellAtlasSection.size() != m_cellEta.size() ) ){
    MSG_FATAL( "FreeRunningHits needs the ATLAS tile cells with their ATLAS identifier fields and the dual readout." );
  }

  return StatusCode::SUCCESS;
}

//!=====================================================================

StatusCode CaloHitMaker::finalize()
{
  return StatusCode::SUCCESS;
}

//!=====================================================================

StatusCode CaloHitMaker::bookHistograms( SG::EventContext &ctx ) const
{
  MSG_DEBUG("Booking histograms...");
  auto store = ctx.getStoreGateSvc();

  store->mkdir(m_histPath);

  // Create the 2D histogram for monitoring purpose
  store->add(new TH2F( "hits_edep", "Hit Energy (Truth); #eta; #phi; Energy [MeV]", m_nEtaBins, m_etaBins.data(), m_nPhiBins, m_phiBins.data() ) );
  
  return StatusCode::SUCCESS;
}

//!=====================================================================

/**
 * @brief Prepare the hit collection before processing steps.
 * 
 * Initializes the CaloHitCollection in the event store. Pre-populates it with 
 * empty hits for every defined readout cell (bin) to be ready for energy accumulation.
 * This ensures that every cell has a corresponding object, even if empty.
 */
StatusCode CaloHitMaker::pre_execute( EventContext &ctx ) const
{
  // One free-running hit list per event, shared by every sampling: the first maker of the sequence creates it.
  if( m_freeRunningHits && !ctx.exist(m_freeRunningListKey) ){
    SG::WriteHandle<CaloFreeRunningHitList> list( m_freeRunningListKey, ctx );
    list.record( std::unique_ptr<CaloFreeRunningHitList>(new CaloFreeRunningHitList()) );
  }

  // Build the CaloHitCollection and attach into the EventContext
  // Create the hit collection into the event context
  SG::WriteHandle<xAOD::CaloHitCollection> collection( m_collectionKey, ctx );
  collection.record( std::unique_ptr<xAOD::CaloHitCollection>(new xAOD::CaloHitCollection()) );

  float deltaEta = std::abs(m_etaBins[1] - m_etaBins[0]);
  float deltaPhi = std::abs(m_phiBins[1] - m_phiBins[0]);

  // Cells given as (r, z) boxes: one hit per cell and phi slice, cells in the order of the table (increasing z),
  // with the eta and delta eta of the table; the local hash follows the same rule as the eta x phi grid.
  if( useCells() ){
    for ( unsigned cell = 0; cell < m_cellEta.size(); ++cell ){
      for ( unsigned phiBin = 0; phiBin < m_nPhiBins; ++phiBin ){
        float phiCenter = m_phiBins[phiBin] + deltaPhi / 2;
        unsigned bin = m_nPhiBins * cell + phiBin;
        auto *hit = new xAOD::CaloHit( m_cellEta[cell], phiCenter, m_cellDeltaEta[cell], deltaPhi, hash(bin),
                                       (CaloSampling)m_sampling,
                                       (Detector)m_detector,
                                       m_bc_duration, m_bcid_start, m_bcid_end );
        hit->setDualReadout( m_tileDualReadout > 0 && isTile() );
        if( !collection->insert( hit->hash(), hit) )
        {
          MSG_FATAL( "It is not possible to include hit hash ("<< hit->hash() << ") into the collection. hash already exist.");
        }
      }
    }
    MSG_DEBUG("Pre_execute done.");
    return StatusCode::SUCCESS;
  }

  //
  // Prepare all sensitive objects like a two dimensional histogram
  //
  for ( unsigned etaBin = 0; etaBin < m_nEtaBins; ++etaBin){ // Rows

    //if (std::abs(m_etaBins[etaBin]) == std::abs(m_etaBins[etaBin+1]))
    //  continue;

    for ( unsigned phiBin = 0; phiBin < m_nPhiBins; ++phiBin){ // Cols

      float etaCenter = m_etaBins[etaBin] + deltaEta / 2;
      float phiCenter = m_phiBins[phiBin] + deltaPhi / 2;
     
      // local hash
      unsigned bin = m_nPhiBins * etaBin + phiBin;
      
      // Create the calorimeter cell
      auto *hit = new xAOD::CaloHit( etaCenter, phiCenter, deltaEta, deltaPhi, hash(bin), 
                                     (CaloSampling)m_sampling,
                                     (Detector)m_detector,
                                     m_bc_duration, m_bcid_start, m_bcid_end );
      if( !collection->insert( hit->hash(), hit) )
      {
        MSG_FATAL( "It is not possible to include hit hash ("<< hit->hash() << ") into the collection. hash already exist.");
      }

    } // Loop over phi bins
  }// Loop over eta bins

  MSG_DEBUG("Pre_execute done.");
  return StatusCode::SUCCESS;
}

//!=====================================================================

/**
 * @brief Execution per Geant4 Step.
 * 
 * Invoked by the framework for every particle step.
 * 1. Checks if the step is within the geometric bounds of this calorimeter volume.
 * 2. Calculates the corresponding eta/phi bin.
 * 3. Retrieves the specific Hit object for that bin.
 * 4. Adds the energy deposit to the Hit.
 */
StatusCode CaloHitMaker::execute( EventContext &ctx , const G4Step *step ) const
{
  SG::ReadHandle<xAOD::CaloHitCollection> collection( m_collectionKey, ctx );

  if( !collection.isValid() ){
    MSG_FATAL("It's not possible to retrieve the CaloHitCollection using this key: " << m_collectionKey);
  }

  // Active medium only: a Geant4 step never crosses a volume boundary, so the material
  // of its pre-step point is the material where the energy was deposited.
  if( m_activeEnergyOnly ){
    const G4Material *material = step->GetPreStepPoint()->GetMaterial();
    if( !material ) return StatusCode::SUCCESS;
    const G4String &materialName = material->GetName();
    if( materialName != "liquidArgon" && materialName != "PLASTIC SCINTILLATOR" )
      return StatusCode::SUCCESS;
  }

  // Get the position
  G4ThreeVector pos = step->GetPreStepPoint()->GetPosition();

  // Apply all necessary transformation (x,y,z) to (eta,phi,r) coordinates
  // Get ATLAS coordinates (in transverse plane xy)
  auto vpos = TVector3( pos.x(), pos.y(), pos.z());
  float eta = vpos.PseudoRapidity();
  float phi = vpos.Phi();
  float radius = vpos.Perp();

  // In plan xy
  if( !(radius > m_rMin && radius <= m_rMax) )
    return StatusCode::SUCCESS;

  if( !(pos.z() > m_zMin && pos.z() <= m_zMax))
    return StatusCode::SUCCESS;

  // Row of the cell: the eta bin of the grid, or the cell of the (r, z) boxes
  int etaBin = useCells() ? findCell(cellRadius(radius, phi), pos.z()) : find(m_etaBins, eta);

  if(etaBin < 0)
    return StatusCode::SUCCESS;

  int phiBin = find(m_phiBins, phi);

  if(phiBin < 0)
    return StatusCode::SUCCESS;

  int bin = m_nPhiBins * etaBin + phiBin;

  xAOD::CaloHit *hit=nullptr;

  if(collection->retrieve(hash(bin), hit)){
    // hit->fill( step );
    const float edep = m_birksLaw ? birksEnergy(step) : (float)step->GetTotalEnergyDeposit();
    if( m_birksLaw )
      hit->fill( step , 1*m_noiseStd, edep );
    else
      hit->fill( step , 1*m_noiseStd); // hit with tof selection sensible by 1*sigma of sampling noise.
    // Dual readout: the same energy shared between the two PMTs by the position across the module
    if( hit->dualReadout() ){
      float w0, w1;
      tilePmtWeights( std::remainder( phi - hit->phi(), 2 * M_PI ), pos.z(), w0, w1 );
      hit->fillPmt( step, edep * w0, edep * w1 );
      // Free-running hits: the same deposit and PMT weights, with the ATLAS identifier and time
      if( m_freeRunningHits && edep > 0 )
        freeRunningTile( ctx, step, etaBin, phiBin, edep, w0, w1 );
    }
  }else{
    MSG_FATAL( "Its not possible to retrieve the hit. Bin ("<< bin << ") not exist");
  }
  
  return StatusCode::SUCCESS;
}

//!=====================================================================

//!=====================================================================

// Weights of the two PMTs of a tile cell for a step at the azimuthal distance phiLocal (rad) from the centre of the
// module, following the ATLAS simulation (Athena, TileGeoG4SDCalc::MakePmtEdepTime):
// - U-shape (TileDualReadout = 1, Ushape = 1 in Athena): the measured response of each PMT across the tile,
//   TileGeoG4SDCalc::Tile_1D_profileRescaled, tables per layer (A, BC, D) for the long barrel and for each extended
//   barrel (TileUShape.h). As in Athena, the PMT 1 reads the table at -phiLocal and the PMT 0 at +phiLocal, and the
//   two weights do not add up exactly to 1.
// - linear (TileDualReadout = 2, Ushape = 0 in Athena): w1 = 0.3 + 0.2 dy1/halfY = 0.5 + 0.2 u, w0 = 1 - w1, with
//   u = tan(phiLocal)/tan(dphi/2) from -1 to 1 across the module (the tile edges taken as radial lines).
// Convention (not checked against the orientation of the ATLAS local axes): the PMT 1 is the one on the side of
// larger phi, on both sides of the detector.
void CaloHitMaker::tilePmtWeights( float phiLocal, float z, float &w0, float &w1 ) const
{
  if( m_tileDualReadout == 1 ){
    // layer: TileCal1/TileExt1 = A, TileCal2/TileExt2 = BC (B), TileCal3/TileExt3 = D
    const int layer = (m_sampling - 5) % 3;
    const int part  = (m_sampling <= 7) ? 0 : ( z > 0 ? 1 : 2 ); // LB, EBA (z > 0), EBC (z < 0)
    const double *t = TileUShape::table( part, layer );
    w1 = TileUShape::amplitude( t, -phiLocal );
    w0 = TileUShape::amplitude( t,  phiLocal );
    return;
  }
  const float halfDeltaPhi = std::abs( m_phiBins[1] - m_phiBins[0] ) / 2;
  float u = std::tan( phiLocal ) / std::tan( halfDeltaPhi );
  u = std::max( -1.f, std::min( 1.f, u ) );
  w1 = 0.5 + 0.2 * u;
  w0 = 1 - w1;
}

//!=====================================================================

// Free-running hits of the ATLAS tile cells, per PMT, with the identifier and the time of the ATLAS simulation.
// Sources (Athena, public, Apache 2.0; copies and report in the F41/F39 notes of the group):
// - identifier: the pmt_id of the ATLAS identifier dictionary (DetectorDescription/IdDictParser/data/
//   IdDictTileCalorimeter.xml), packed in 64 bits as the IdDict code does (IdDictDictionary, IdentifierField):
//   subdet (Tile = index 2) in bits 63-61, section 60-58, side 57-54 (index: -1 -> 0, +1 -> 2), module 53-46,
//   tower 45-40, sampling 39-36, pmt 35-34 (adc 33-32 = 0); checked against Calorimeter/CaloIdentifier/share/TileID_test.ref.
//   The PMT on the side of larger phi ("up") is pmt 1, the other pmt 0 (the convention of tilePmtWeights); D0 is side +1.
// - module m covers phi in [m, m+1] * 2 pi / 64 (TileDetDescr); side +1 for z > 0.
// - time (TileGeoG4SDCalc::MakePmtEdepTime, DoTOFCorrection = true, refraction index 1.59 of tile and fibre): the global
//   time at the post-step point, plus [r_cell (n - 2 / sin(theta_cell)) + |x| (1 - n sin(theta_hit))] / c, where r_cell and
//   theta_cell are those of the centre of the cell and x the post-step position; each PMT adds
//   n (dy - r_cell tan(pi/64)) / c, dy the distance of the step to the edge of the tile on the side of the other PMT
//   (here: the edges of the module at +-pi/64 around its centre, from the pre-step position, as for the PMT weights).
//   If the time of a PMT passes TimeCut = 350.5 ns, both go to the late hit time 99995 ns without the PMT delays; the
//   energy is kept.
// - time binning (option FreeRunningBinning; TileGeoG4SDCalc::deltaT, DeltaTHit = {0.5, -75.25, 75.25, 5.}): bins of
//   0.5 ns for -75.25 < t < 75.25 ns, 5 ns elsewhere (see CaloFreeRunningHitList::add).
void CaloHitMaker::freeRunningTile( SG::EventContext &ctx, const G4Step *step, int cell, int phiBin, float edep,
                                    float w0, float w1 ) const
{
  SG::ReadHandle<CaloFreeRunningHitList> list( m_freeRunningListKey, ctx );
  if( !list.isValid() ){
    MSG_FATAL("It's not possible to retrieve the CaloFreeRunningHitList using this key: " << m_freeRunningListKey);
  }
  const double n = 1.59, timeCut = 350.5, lateHitTime = 99995., tanPi64 = std::tan(M_PI / 64);
  const double cLight = CLHEP::c_light;   // mm/ns

  // identifier fields
  const double twoPi = 2 * M_PI, modWidth = twoPi / 64;
  double phiCentre = 0.5 * ( m_phiBins[phiBin] + m_phiBins[phiBin+1] );
  double phi0 = phiCentre < 0 ? phiCentre + twoPi : phiCentre;
  int module = (int)std::floor( phi0 / modWidth ) % 64;
  double moduleCentre = std::remainder( (module + 0.5) * modWidth, twoPi );
  int section  = m_cellAtlasSection[cell];
  int tower    = m_cellAtlasTower[cell];
  int sampling = m_cellAtlasSampling[cell];
  int side     = ( m_cellZCentre[cell] >= 0 || (section == 1 && sampling == 2 && tower == 0) ) ? +1 : -1;
  auto pmtId = [&]( int pmt ) -> long long {
    unsigned long long id = (2ULL << 61) | ((unsigned long long)section << 58) | ((unsigned long long)(side > 0 ? 2 : 0) << 54)
                          | ((unsigned long long)module << 46) | ((unsigned long long)tower << 40)
                          | ((unsigned long long)sampling << 36) | ((unsigned long long)pmt << 34);
    return (long long)id;
  };

  // time of the deposit (ns)
  const G4ThreeVector post = step->GetPostStepPoint()->GetPosition();
  const G4ThreeVector pre  = step->GetPreStepPoint()->GetPosition();
  double tglobal = step->GetPostStepPoint()->GetGlobalTime() / ns;
  double rCell = m_cellRCentre[cell], zCell = m_cellZCentre[cell];
  double sinThCell = rCell / std::sqrt( rCell*rCell + zCell*zCell );
  double cosThHit = post.cosTheta();
  double t = tglobal + ( rCell * ( n - 2.0 / sinThCell ) + post.mag() * ( 1.0 - n * std::sqrt( 1.0 - cosThHit*cosThHit ) ) ) / cLight;
  // PMT delays: light across the tile to each edge; local frame of the module from the pre-step position
  double phiLocal = std::remainder( pre.phi() - moduleCentre, twoPi );
  double rPre = pre.perp();
  double xLocal = rPre * std::cos(phiLocal), yLocal = rPre * std::sin(phiLocal);
  double halfWidth = xLocal * tanPi64;
  double dyUp = halfWidth - yLocal;     // to the edge on the side of larger phi (PMT 1 reads there)
  double dyDown = halfWidth + yLocal;
  double dtUp   = n * ( dyUp   - rCell * tanPi64 ) / cLight;
  double dtDown = n * ( dyDown - rCell * tanPi64 ) / cLight;
  if( t + dtUp > timeCut || t + dtDown > timeCut ){
    t = lateHitTime; dtUp = dtDown = 0.0;
  }
  double delta = m_freeRunningBinning ? ( ( t > -75.25 && t < 75.25 ) ? 0.5 : 5.0 ) : 0.0;

  double eta = m_cellEta[cell];
  if( w1 * edep > 0 ) list->add( pmtId(1), w1 * edep, t + dtUp,   tglobal, eta, moduleCentre, sampling, side, module, tower, delta );
  if( w0 * edep > 0 ) list->add( pmtId(0), w0 * edep, t + dtDown, tglobal, eta, moduleCentre, sampling, side, module, tower, delta );
}

//!=====================================================================

// Energy of the step after Birks' law, following the ATLAS simulation (Athena):
// - plastic scintillator (TileGeoG4SDCalc::BirkLaw): E / (1 + kB dE/dx), kB = 0.02002 g/(MeV cm2)
//   ("value updated for G4 10.6.p03"), second-order term 0, only for charged particles,
//   kB scaled by 7.2/12.6 for charge above 1;
// - liquid argon (LArG4BirksLaw): E (1 + k/F 1.51) / (1 + k/F (dE/dx)/rho), k = 0.05832,
//   rho = 1.396 g/cm3, F = 10 kV/cm (the constant field of the ATLAS barrel, presamplers and HEC),
//   with the correction for heavily ionising particles above 969 MeV/cm.
// Other materials are left unchanged.
float CaloHitMaker::birksEnergy( const G4Step *step ) const
{
  const double edep = step->GetTotalEnergyDeposit();
  const double length = step->GetStepLength();
  const G4Material *material = step->GetPreStepPoint()->GetMaterial();
  if( !material || edep <= 0 || length <= 0 ) return (float)edep;
  const G4String &name = material->GetName();

  if( name == "PLASTIC SCINTILLATOR" ){
    const double charge = step->GetPreStepPoint()->GetCharge();
    if( charge == 0. ) return (float)edep;
    double kB = 0.02002 * g / (MeV * cm2);
    if( std::abs(charge) > 1.0 ) kB *= 7.2 / 12.6;
    const double dedx = edep / length / material->GetDensity();
    return (float)( edep / (1. + kB * dedx) );
  }

  if( name == "liquidArgon" ){
    const double dE = edep / MeV;
    const double dX = length / cm;
    if( dX < 1e-5 ) return (float)edep;
    const double kOverField = 0.05832 / 10.;
    double dEdX = dE / dX;
    const double corrected = dE * (1 + kOverField * 1.51) / (1 + kOverField * dEdX / 1.396);
    if( dEdX > 12000.0 ) dEdX = 12000.0;
    const double kHIP = (dEdX > 969.) ? 0.000754 * dEdX + 0.2692 : 1.0;
    return (float)( corrected * kHIP * MeV );
  }

  return (float)edep;
}

//!=====================================================================

// standlone execute
StatusCode CaloHitMaker::execute( EventContext &/*ctx*/, int /*evt*/ ) const
{
  MSG_ERROR("This method can not be execute in standalone mode.");
  return StatusCode::FAILURE;
}

//!=====================================================================

StatusCode CaloHitMaker::post_execute( EventContext &/*ctx*/ ) const
{
  return StatusCode::SUCCESS;
}

//!=====================================================================

StatusCode CaloHitMaker::fillHistograms( EventContext &ctx ) const
{
  auto store = ctx.getStoreGateSvc();
  SG::ReadHandle<xAOD::CaloHitCollection> collection( m_collectionKey, ctx );
 
  if( !collection.isValid() ){
    MSG_FATAL("It's not possible to retrieve the CaloHitCollection using this key: " << m_collectionKey);
  }

  store->cd(m_histPath);

  
  for ( const auto& p : **collection.ptr() ){ 

    const auto *hit = p.second;
    
    {// Fill truth energy 2D histograms
      int x = store->hist2("hits_edep")->GetXaxis()->FindBin(hit->eta());
      int y = store->hist2("hits_edep")->GetYaxis()->FindBin(hit->phi());
      int bin = store->hist2("hits_edep")->GetBin(x,y,0);
      float energy = store->hist2("hits_edep")->GetBinContent( bin );
      store->hist2("hits_edep")->SetBinContent( bin, (energy + hit->edep()) );
    }
  }

  return StatusCode::SUCCESS;
}

//!=====================================================================

// Cell of the (r, z) boxes that contains the point (rmin < r <= rmax, zmin < z <= zmax), or -1
int CaloHitMaker::findCell( float radius, float z ) const
{
  for ( unsigned box = 0; box < m_cellBoxIndex.size(); ++box ){
    if( radius > m_cellBoxRMin[box] && radius <= m_cellBoxRMax[box] &&
        z > m_cellBoxZMin[box] && z <= m_cellBoxZMax[box] )
      return m_cellBoxIndex[box];
  }
  return -1;
}

//!=====================================================================

// Radius used to find the (r, z) box of a step. With CellRadiusModules = N > 0: the distance to the beam line measured
// on the axis of the module that contains phi, r cos(phi - phi_c), with modules of 2 pi / N from phi = 0 and phi_c
// the centre of the module. This is moduleY of the ATLAS HEC simulation (Athena, LArCalorimeter/LArG4/LArG4HEC/src/
// HECGeometry.cc, CalculateIdentifier: |y| of the step in the frame of its module), for 32 modules with the first
// one starting at phi = 0 (phi0 = 0 of the HEC in the ATLAS identifier dictionary). N = 0 (default): r.
float CaloHitMaker::cellRadius( float radius, float phi ) const
{
  if( m_cellRadiusModules <= 0 ) return radius;
  const double twoPi = 2 * M_PI, width = twoPi / m_cellRadiusModules;
  double phi0 = phi < 0 ? phi + twoPi : phi;
  int module = (int)std::floor( phi0 / width );
  if( module >= m_cellRadiusModules ) module = m_cellRadiusModules - 1;
  return (float)( radius * std::cos( phi0 - (module + 0.5) * width ) );
}

//!=====================================================================

int CaloHitMaker::find( const std::vector<float> &vec, float value) const
{
  auto binIterator = std::adjacent_find( vec.begin(), vec.end(), [=](float left, float right){ return left < value and value <= right; }  );
  if ( binIterator == vec.end() ) return -1;
  return  binIterator - vec.begin();
}

//!=====================================================================
/*
class CaloSampling(EnumStringification):
    PSB       = 0
    PSE       = 1
    EMB1      = 2
    EMB2      = 3
    EMB3      = 4
    TileCal1  = 5
    TileCal2  = 6
    TileCal3  = 7
    TileExt1  = 8
    TileExt2  = 9
    TileExt3  = 10
    EMEC1     = 11
    EMEC2     = 12
    EMEC3     = 13
    HEC1      = 14
    HEC2      = 15
    HEC3      = 16
    HEC0      = 60   (ATLAS-like HEC, --atlas-hec)
*/
unsigned long int CaloHitMaker::hash(unsigned bin) const
{
  // NOTE: Max of 4,294,967,295 for unsigned long int (10 digits)

  // two for sammpling+side (0-99), one for segment (0-9) and six for bin index (0-999999), 9 digits
  if(m_zMin > 0){ // Right side 
    return ( (m_sampling + 1*17 ) * 1e7 + m_segment * 1e6 + bin);
  }else if(m_zMax < 0){// Left side
    return ( (m_sampling + 2*17 ) * 1e7 + m_segment * 1e6 + bin);
  }else{ // Center (barrel)
    return ( (m_sampling + 0*17 ) * 1e7 + m_segment * 1e6 + bin);
  }
}

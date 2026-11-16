
#include "HardFilter.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include <cmath>
#include <vector>

using namespace generator;

namespace {

  // Thresholds of DSID 800831 (GeV)
  const double kJetR        = 0.6;
  const double kJetPtMin    = 35.;
  const double kPhotonPtMin = 8.;
  const double kBPtMin      = 5.;
  const double kLeptonPtMin = 4.;
  const double kEtaMax      = 4.5;

  bool isStable( const HepMC3::ConstGenParticlePtr &p ){ return p->status() == 1; }

  bool isNeutrino( int apid ){ return apid == 12 || apid == 14 || apid == 16 || apid == 18; }

  // MC::isWeaklyDecayingBHadron (TruthUtils/AtlasPID.h)
  bool isWeaklyDecayingBHadron( int pid ){
    const int a = std::abs(pid);
    return a == 511  || a == 521  || a == 531  || a == 541  ||
           a == 5122 || a == 5132 || a == 5232 || a == 5332 ||
           a == 5142 || a == 5242 || a == 5342 || a == 5442 ||
           a == 5512 || a == 5522 || a == 5532 || a == 5542;
  }

  // DirectPhotonFilter: a photon is not prompt when any parent has |pdg id| > 100
  bool fromHadron( const HepMC3::ConstGenParticlePtr &p ){
    const auto vtx = p->production_vertex();
    if ( !vtx ) return false;
    for ( const auto &parent : vtx->particles_in() ){
      if ( std::abs(parent->pid()) > 100 ) return true;
    }
    return false;
  }

}


int HardFilter::evaluate( const HepMC3::GenEvent &evt, double *leadingJetPt )
{
  int pass = 0;
  bool hasB = false, hasLepton = false;
  std::vector<fastjet::PseudoJet> inputs;

  for ( const auto &p : evt.particles() ){
    const int pid  = p->pid();
    const int apid = std::abs(pid);
    const auto &mom = p->momentum();
    const double pt = mom.perp();

    // B hadrons are checked at any status (they decay)
    if ( !hasB && isWeaklyDecayingBHadron(pid) && pt > kBPtMin && std::abs(mom.eta()) < kEtaMax )
      hasB = true;

    if ( !isStable(p) ) continue;

    // Leptons: leading stable electron or muon
    if ( !hasLepton && (apid == 11 || apid == 13) && pt > kLeptonPtMin && std::abs(mom.eta()) <= kEtaMax )
      hasLepton = true;

    // Prompt photons
    if ( !(pass & Photon) && pid == 22 && pt > kPhotonPtMin && std::abs(mom.eta()) <= kEtaMax && !fromHadron(p) )
      pass |= Photon;

    // Truth jet inputs: stable, no neutrinos, no muons, no zero energy photons
    if ( isNeutrino(apid) || apid == 13 ) continue;
    if ( pid == 22 && mom.e() <= 0 ) continue;
    inputs.emplace_back( mom.px(), mom.py(), mom.pz(), mom.e() );
  }

  if ( hasB && hasLepton ) pass |= BLepton;

  double lead = 0;
  if ( !inputs.empty() ){
    fastjet::JetDefinition def( fastjet::antikt_algorithm, kJetR );
    fastjet::ClusterSequence cs( inputs, def );
    for ( const auto &j : cs.inclusive_jets( 0. ) ) lead = std::max( lead, j.pt() );
  }
  if ( lead > kJetPtMin ) pass |= Jet;
  if ( leadingJetPt ) *leadingJetPt = lead;

  return pass;
}

#ifndef Evtgen_HardFilter_h
#define Evtgen_HardFilter_h

#include "HepMC3/GenEvent.h"

namespace generator{

  /**
   * @brief Selection of the ATLAS high-pT minimum bias sample used for pileup
   * (mc.Py8EG_minbias_inelastic_highjetphotonlepton, DSID 800831, in MCJobOptions/800xxx).
   *
   * An inelastic event is "hard" when at least one of the following holds:
   *   - an anti-kt R = 0.6 truth jet with pT > 35 GeV (QCDTruthJetFilter, no eta cut);
   *   - a stable photon with pT > 8 GeV and |eta| < 4.5 that does not come from a hadron decay
   *     (DirectPhotonFilter);
   *   - a weakly decaying B hadron with pT > 5 GeV and |eta| < 4.5 (HeavyFlavorHadronFilter)
   *     and a stable electron or muon with pT > 4 GeV and |eta| < 4.5 (LeptonFilter).
   * The truth jets are built as the ATLAS generator jets (CopyTruthJetParticles defaults):
   * stable particles, without neutrinos and muons.
   * The ATLAS sample decays the B hadrons with EvtGen; here Pythia decays them. The B hadron
   * selection does not depend on the decay, since the B hadron is in the record before it decays.
   */
  namespace HardFilter{

    // Bit pattern of the criteria that hold (0 when the event is not hard)
    enum Criterion { Jet = 1, Photon = 2, BLepton = 4 };

    int evaluate( const HepMC3::GenEvent &evt, double *leadingJetPt=nullptr );

  }
}

#endif

#ifndef Evtgen_Pileup_h
#define Evtgen_Pileup_h

#include "GenKernel/IAlgorithm.h"
#include "GenKernel/IGenerator.h"
#include <memory>
#include <string>

namespace HepMC3 { class WriterRootTree; }

namespace generator{

  class Pileup: public IAlgorithm 
  {
    public:
  
      Pileup(const std::string , IGenerator *);
      ~Pileup();
  
      virtual StatusCode initialize() override;
      virtual StatusCode execute( Event& ) override;
      virtual StatusCode finalize() override;
  
    private:
      
      int m_select;
      int m_bc_id_start;
      int m_bc_id_end;
      float m_etaMax;
      float m_nPileupAvg;
      float m_nPileupSigma;
      float m_nPileupPerBunch;
      float m_delta_eta;
      float m_delta_phi;
      float m_ptMinCharged;
      float m_ptMinNeutral;
      int m_hardFilter;          // 0: off, 1: mark the hard events (DSID 800831 selection), 2: keep only hard events
      std::string m_truthFile;   // HepMC3 ROOT file with the full generator record of every event written (empty: none)
      std::shared_ptr<HepMC3::WriterRootTree> m_truthWriter; //!
      long m_nGenerated = 0; //!
      long m_nKept = 0; //!

  };
}
#endif

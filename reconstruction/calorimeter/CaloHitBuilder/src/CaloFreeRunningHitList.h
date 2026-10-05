#ifndef CaloFreeRunningHitList_h
#define CaloFreeRunningHitList_h

#include "GaugiKernel/DataHandle.h"
#include <cmath>
#include <map>
#include <vector>

/**
 * @class CaloFreeRunningHitList
 * @brief Per-event list of the tile hits in the layout of the ATLAS HITS ntuple (one entry per PMT and time bin, or per
 *        PMT and Geant4 step without binning), filled by CaloHitMaker when simu_trf.py runs with --free-running-hits and
 *        written by CaloFreeRunningHitWriter.
 *
 * With binning (bin width delta > 0), add() follows TileSimHit::add of the ATLAS simulation (Athena,
 * TileCalorimeter/TileSimEvent/src/TileSimHit.cxx): a deposit of a PMT goes to an existing bin of the same PMT when
 * |t - t_bin| < delta/2, otherwise it opens a new bin at t = delta * nearbyint(t/delta). The energies of a bin are summed;
 * the extra tglobal of a bin is the energy weighted mean of the Geant4 global times of its deposits.
 *
 * The event context only gives const access to recorded objects, so the vectors are mutable and add() is const (the
 * same situation as CaloHitMaker filling the CaloHit objects of a const collection). Each event is simulated by a single
 * thread, so the list of an event is never filled concurrently.
 */
class CaloFreeRunningHitList : public SG::DataHandle
{
  public:

    CaloFreeRunningHitList()=default;
    ~CaloFreeRunningHitList()=default;

    void add( long long cellID, double energy, double time, double tglobal, double eta, double phi,
              int sampling, int side, int module, int tower, double delta ) const
    {
      if( delta > 0 ){
        auto &bins = m_bins[cellID];
        for( size_t i : bins ){
          if( std::fabs(time - m_time[i]) < delta/2. ){
            m_energy[i]    += energy;
            m_sumETglobal[i] += energy * tglobal;
            return;
          }
        }
        bins.push_back( m_cellID.size() );
        push( cellID, energy, delta * std::nearbyint(time/delta), tglobal, eta, phi, sampling, side, module, tower );
      }else{
        push( cellID, energy, time, tglobal, eta, phi, sampling, side, module, tower );
      }
    }

    size_t size() const { return m_cellID.size(); }

    const std::vector<long long>& cellID()   const { return m_cellID;   }
    const std::vector<double>&    energy()   const { return m_energy;   }
    const std::vector<double>&    time()     const { return m_time;     }
    const std::vector<double>&    eta()      const { return m_eta;      }
    const std::vector<double>&    phi()      const { return m_phi;      }
    const std::vector<int>&       sampling() const { return m_sampling; }
    const std::vector<int>&       side()     const { return m_side;     }
    const std::vector<int>&       module()   const { return m_module;   }
    const std::vector<int>&       tower()    const { return m_tower;    }
    /*! energy weighted mean of the Geant4 global times of the deposits of each entry [ns] */
    std::vector<double> tglobal() const
    {
      std::vector<double> t( m_cellID.size() );
      for( size_t i=0; i<t.size(); ++i ) t[i] = m_energy[i] > 0 ? m_sumETglobal[i] / m_energy[i] : m_firstTglobal[i];
      return t;
    }

  private:

    void push( long long cellID, double energy, double time, double tglobal, double eta, double phi,
               int sampling, int side, int module, int tower ) const
    {
      m_cellID.push_back(cellID);   m_energy.push_back(energy);   m_time.push_back(time);
      m_sumETglobal.push_back(energy * tglobal);   m_firstTglobal.push_back(tglobal);
      m_eta.push_back(eta);   m_phi.push_back(phi);
      m_sampling.push_back(sampling);   m_side.push_back(side);   m_module.push_back(module);   m_tower.push_back(tower);
    }

    /*! ATLAS PMT identifier (64-bit compact Identifier of the ATLAS identifier dictionary) */
    mutable std::vector<long long> m_cellID;
    /*! visible energy [MeV] */
    mutable std::vector<double>    m_energy;
    /*! hit time of the ATLAS simulation (TileGeoG4SDCalc) [ns]; bin centre with binning */
    mutable std::vector<double>    m_time;
    mutable std::vector<double>    m_sumETglobal;
    mutable std::vector<double>    m_firstTglobal;
    /*! centre of the cell in eta, centre of the module in phi */
    mutable std::vector<double>    m_eta;
    mutable std::vector<double>    m_phi;
    /*! ATLAS fields: sampling (0 A, 1 BC/B, 2 D), side (+1 for z > 0), module (0-63), tower */
    mutable std::vector<int>       m_sampling;
    mutable std::vector<int>       m_side;
    mutable std::vector<int>       m_module;
    mutable std::vector<int>       m_tower;
    /*! entries (bins) of each PMT, for the binning */
    mutable std::map<long long, std::vector<size_t>> m_bins;
};

#endif

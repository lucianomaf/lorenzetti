"""The C++ routines of the generator, declared once in ROOT: the reader of entry ranges of a hits tree of the library, the AR(4) recursion
of the LAr noise, the exponential sum of the Tile tail and the writer of the hits tree of the output. The recursions keep the order of the
operations of a plain loop and do not fuse multiplications and additions (fp contract off): the same numbers as the loop in Python."""

_ROOT = []

CODE = r"""
#include <TFile.h>
#include <TTree.h>
#include <memory>
#include <vector>

struct LztSlices { std::vector<int> ev; std::vector<float> e, t; };

// the entries [lo[i], hi[i]) of the tree hits (event int, energy float, time float), concatenated in the order of the ranges
LztSlices lzt_read_ranges(const char *fn, const std::vector<Long64_t> &lo, const std::vector<Long64_t> &hi) {
  LztSlices s; Long64_t n = 0;
  for (size_t i = 0; i < lo.size(); ++i) n += hi[i] - lo[i];
  s.ev.resize(n); s.e.resize(n); s.t.resize(n);
  std::unique_ptr<TFile> f(TFile::Open(fn));
  TTree *t = f->Get<TTree>("hits");
  int ev; float e, tm;
  t->SetBranchAddress("event", &ev); t->SetBranchAddress("energy", &e); t->SetBranchAddress("time", &tm);
  Long64_t k = 0;
  for (size_t i = 0; i < lo.size(); ++i)
    for (Long64_t j = lo[i]; j < hi[i]; ++j) { t->GetEntry(j); s.ev[k] = ev; s.e[k] = e; s.t[k] = tm; ++k; }
  t->ResetBranchAddresses();
  return s;
}

// y[i] = a1 y[i-1] + a2 y[i-2] + a3 y[i-3] + a4 y[i-4] + e[i] from i = 4, y[0..3] = 0
void lzt_ar4(const double *e, Long64_t n, double a1, double a2, double a3, double a4, double *y) {
#pragma clang fp contract(off)
  for (Long64_t i = 0; i < n && i < 4; ++i) y[i] = 0.0;
  for (Long64_t i = 4; i < n; ++i) y[i] = a1 * y[i - 1] + a2 * y[i - 2] + a3 * y[i - 3] + a4 * y[i - 4] + e[i];
}

// acc = acc f + inj[i]; y[i] += acc
void lzt_exp_sum(const double *inj, Long64_t n, double f, double *y) {
#pragma clang fp contract(off)
  double acc = 0.0;
  for (Long64_t i = 0; i < n; ++i) { acc = acc * f + inj[i]; y[i] += acc; }
}

// the hits tree of the output, in the current directory
void lzt_write_hits(Long64_t n, const int *ch, const int *smp, const float *e, const float *tm, const int *bc, const int *ev, const int *cp) {
  TTree *t = new TTree("hits", "the hits of the kept samples");
  int a, b, c, d, g; float x, y;
  t->Branch("channel", &a, "channel/I"); t->Branch("sample", &b, "sample/I"); t->Branch("visible_energy_mev", &x, "visible_energy_mev/F");
  t->Branch("time_ns", &y, "time_ns/F"); t->Branch("origin_bc", &c, "origin_bc/I"); t->Branch("event", &d, "event/I");
  t->Branch("copy", &g, "copy/I");
  for (Long64_t i = 0; i < n; ++i) { a = ch[i]; b = smp[i]; x = e[i]; y = tm[i]; c = bc[i]; d = ev[i]; g = cp[i]; t->Fill(); }
  t->Write();
  delete t;
}
"""


def load():
    """ROOT with the routines declared"""
    if not _ROOT:
        import ROOT
        ROOT.gInterpreter.Declare(CODE)
        _ROOT.append(ROOT)
    return _ROOT[0]

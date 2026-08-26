// ---------------------------------------------------------------------------
// hist_merge.cxx — combine the per-shard histogram files written by the
// *_hist stages into one file the *_plot stages can read unchanged.
//
//   ./hist_merge <output.root> <input1.root> <input2.root> ...
//
// WHY NOT hadd. `hadd` merges anything with a Merge() method and, for
// everything else, keeps the FIRST occurrence. Every scalar HistWriter emits
// is a TParameter -- event counts, rejection tallies, the pT accumulators
// rpt_v5_plot turns into untimed-fraction percentages -- and TParameter has no
// Merge(). hadd would therefore produce a file whose histograms are the sum of
// N shards but whose event counts are one shard's, silently and with no error.
// Every one of those scalars is an additive accumulator (a count, or a sum of
// pT that the plot stage only ever divides by another summed scalar), so the
// correct merge is to add them -- which is what this does.
//
// Types handled:
//   TH1 and everything deriving from it (TH1D/TH2D/TProfile/TProfile2D) -> Add
//   TParameter<Long64_t> / TParameter<Double_t>                          -> sum
//   TObjString (meta_energy_label)                                       -> must
//        agree across inputs; kept once. A disagreement means the shards came
//        from different samples, which is a merge that should fail loudly.
//
// Correctness of sharding itself: shards partition the SORTED file list
// round-robin (see setupChain), so every event is processed exactly once
// across the set. All histogram fills in this codebase are unweighted, so bin
// contents are integer-valued sums and the merged result is bit-identical to
// an unsharded run -- not merely close. util/scratch/compare_hists.py checks
// that, and CLAUDE.md records the validation.
// ---------------------------------------------------------------------------

#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TObjString.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TError.h>

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace {

// Everything read out of a shard file, keyed by name, detached from any TFile
// so it survives that file being closed.
struct Merged {
  std::map<std::string, TH1*>                    hists;
  std::map<std::string, Long64_t>                longs;
  std::map<std::string, Double_t>                doubles;
  std::map<std::string, std::string>             strings;
  std::vector<std::string>                       order;   // first-seen order
  std::map<std::string, char>                    kind;    // 'h','l','d','s'

  void note(const std::string& n, char k) {
    if (!kind.count(n)) { order.push_back(n); kind[n] = k; }
  }
};

bool absorb(const std::string& path, Merged& m, bool first) {
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
  if (!f || f->IsZombie()) {
    std::cerr << "hist_merge: cannot open '" << path << "'\n";
    return false;
  }

  TIter next(f->GetListOfKeys());
  while (auto* key = static_cast<TKey*>(next())) {
    const std::string name = key->GetName();
    std::unique_ptr<TObject> obj(key->ReadObj());
    if (!obj) continue;

    if (auto* h = dynamic_cast<TH1*>(obj.get())) {
      m.note(name, 'h');
      if (first || !m.hists.count(name)) {
        h->SetDirectory(nullptr);          // detach: outlives this TFile
        m.hists[name] = static_cast<TH1*>(obj.release());
      } else {
        // TH1::Add is correct for TProfile/TProfile2D too -- they override it
        // to combine sums and sums-of-squares rather than bin means.
        if (!m.hists[name]->Add(h)) {
          std::cerr << "hist_merge: Add failed for '" << name << "' in " << path
                    << " (binning mismatch?)\n";
          return false;
        }
      }
    } else if (auto* pl = dynamic_cast<TParameter<Long64_t>*>(obj.get())) {
      m.note(name, 'l');
      m.longs[name] += pl->GetVal();
    } else if (auto* pd = dynamic_cast<TParameter<Double_t>*>(obj.get())) {
      m.note(name, 'd');
      m.doubles[name] += pd->GetVal();
    } else if (auto* os = dynamic_cast<TObjString*>(obj.get())) {
      m.note(name, 's');
      const std::string v = os->GetString().Data();
      auto it = m.strings.find(name);
      if (it == m.strings.end()) m.strings[name] = v;
      else if (it->second != v) {
        std::cerr << "hist_merge: '" << name << "' differs between inputs:\n"
                  << "  have: " << it->second << "\n  " << path << ": " << v << "\n"
                  << "Refusing to merge shards that disagree on run metadata.\n";
        return false;
      }
    } else {
      std::cerr << "hist_merge: unhandled type " << obj->ClassName()
                << " for key '" << name << "' -- refusing to drop it silently.\n";
      return false;
    }
  }
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "usage: hist_merge <output.root> <input1.root> [input2.root ...]\n";
    return 1;
  }
  // Merged histograms are held detached from any file; without this ROOT would
  // also try to attach freshly-read ones to whichever file is current.
  TH1::AddDirectory(kFALSE);
  gErrorIgnoreLevel = kWarning;

  const std::string outPath = argv[1];
  std::vector<std::string> inputs(argv + 2, argv + argc);

  Merged m;
  for (size_t i = 0; i < inputs.size(); ++i) {
    if (!absorb(inputs[i], m, i == 0)) return 1;
    std::cout << "  + " << inputs[i] << '\n';
  }

  std::unique_ptr<TFile> out(TFile::Open(outPath.c_str(), "RECREATE"));
  if (!out || out->IsZombie()) {
    std::cerr << "hist_merge: cannot create '" << outPath << "'\n";
    return 1;
  }
  out->cd();
  size_t nH = 0, nS = 0;
  for (const auto& name : m.order) {
    switch (m.kind[name]) {
      case 'h': m.hists[name]->Write(name.c_str());                          ++nH; break;
      case 'l': TParameter<Long64_t>(name.c_str(), m.longs[name]).Write();   ++nS; break;
      case 'd': TParameter<Double_t>(name.c_str(), m.doubles[name]).Write(); ++nS; break;
      case 's': TObjString(m.strings[name].c_str()).Write(name.c_str());     ++nS; break;
    }
  }
  out->Close();

  std::cout << "hist_merge: wrote " << outPath << "  ("
            << nH << " histograms, " << nS << " scalars, from "
            << inputs.size() << " inputs)\n";
  return 0;
}

// Check if jets are pT-sorted in the ROOT files
// Run with: root -l -q -b check_jet_sorting.cxx

R__ADD_INCLUDE_PATH(/opt/homebrew/opt/boost/include)
R__ADD_LIBRARY_PATH(/opt/homebrew/opt/boost/lib)
R__LOAD_LIBRARY(libboost_filesystem)

#include "src/clustering_functions.h"
#include "src/event_processing.h"
#include <iostream>

using namespace MyUtl;

void check_jet_sorting(const std::string& fileNumber, int maxEvents) {

  // Setup input chain
  TChain chain("ntuple");
  setupChain(chain, fileNumber);

  if (chain.GetEntries() == 0) {
    std::cerr << "No entries found in chain!" << std::endl;
    return;
  }

  // Setup reader and branches
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  int totalEvents = 0;
  int eventsSorted = 0;
  int eventsUnsorted = 0;
  int eventsWithMultipleJets = 0;

  std::cout << "\n========================================" << std::endl;
  std::cout << "Jet pT Sorting Check" << std::endl;
  std::cout << "========================================\n" << std::endl;

  // Event loop
  while (reader.Next() && totalEvents < maxEvents) {
    totalEvents++;

    int nJets = branch.topoJetPt.GetSize();

    if (nJets < 2) continue;

    eventsWithMultipleJets++;

    bool isSorted = true;

    // Print first few events in detail
    if (eventsWithMultipleJets <= 5) {
      std::cout << "Event " << totalEvents << " (nJets = " << nJets << "):" << std::endl;
      for (int i = 0; i < nJets && i < 6; ++i) {
        std::cout << "  Jet " << i << ": pT = " << branch.topoJetPt[i]
                  << " GeV, eta = " << branch.topoJetEta[i] << std::endl;
      }
    }

    // Check if jets are pT-sorted (descending order)
    for (int i = 0; i < nJets - 1; ++i) {
      if (branch.topoJetPt[i] < branch.topoJetPt[i+1]) {
        isSorted = false;
        if (eventsUnsorted < 3) {
          std::cout << "\n!!! UNSORTED EVENT " << totalEvents << " !!!" << std::endl;
          std::cout << "  Jet " << i << ": pT = " << branch.topoJetPt[i] << " GeV" << std::endl;
          std::cout << "  Jet " << i+1 << ": pT = " << branch.topoJetPt[i+1] << " GeV" << std::endl;
          std::cout << "  (Jet " << i+1 << " has higher pT than Jet " << i << ")" << std::endl;
        }
        break;
      }
    }

    if (isSorted) {
      eventsSorted++;
    } else {
      eventsUnsorted++;
    }

    if (eventsWithMultipleJets <= 5) {
      std::cout << "  Sorted: " << (isSorted ? "YES" : "NO") << "\n" << std::endl;
    }
  }

  // Print summary
  std::cout << "\n========================================" << std::endl;
  std::cout << "Jet Sorting Summary" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "Total events processed: " << totalEvents << std::endl;
  std::cout << "Events with 2+ jets: " << eventsWithMultipleJets << std::endl;
  std::cout << "Events with sorted jets: " << eventsSorted
            << " (" << 100.0*eventsSorted/eventsWithMultipleJets << "%)" << std::endl;
  std::cout << "Events with unsorted jets: " << eventsUnsorted
            << " (" << 100.0*eventsUnsorted/eventsWithMultipleJets << "%)" << std::endl;

  if (eventsUnsorted == 0) {
    std::cout << "\n✓ JETS ARE pT-SORTED (descending order)" << std::endl;
    std::cout << "  Leading jet = topoJetPt[0]" << std::endl;
    std::cout << "  Subleading jet = topoJetPt[1]" << std::endl;
  } else {
    std::cout << "\n✗ JETS ARE NOT pT-SORTED" << std::endl;
    std::cout << "  Need to implement sorting!" << std::endl;
  }
  std::cout << "========================================\n" << std::endl;
}

void check_jet_sorting() {
  check_jet_sorting("000001", 1000);
}

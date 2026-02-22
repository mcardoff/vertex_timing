// Test VBF Signal Region Selection Functions
// Demonstrates usage of the new VBF selection methods

#include "src/clustering_functions.h"
#include "src/event_processing.h"
#include <iostream>

using namespace MyUtl;

void test_vbf_selection(const std::string& fileNumber = "000001", int maxEvents = 100) {

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

  // Counters for VBF selection
  int totalEvents = 0;
  int passMultiplicity = 0;
  int passLeadPt = 0;
  int passDeltaPhi = 0;
  int passOppositeHemi = 0;
  int passDeltaEta = 0;
  int passMjj = 0;
  int passPtSum = 0;
  int passForward = 0;
  int passFSR = 0;
  int passFullSelection = 0;

  std::cout << "\n========================================" << std::endl;
  std::cout << "VBF H->Invisible Selection Test" << std::endl;
  std::cout << "========================================\n" << std::endl;

  // Event loop
  while (reader.Next() && totalEvents < maxEvents) {
    totalEvents++;

    // Apply individual cuts and count
    if (branch.passJetMultiplicity()) passMultiplicity++;
    if (branch.passLeadingJetPt()) passLeadPt++;
    if (branch.passJetDeltaPhi()) passDeltaPhi++;
    if (branch.passOppositeHemispheres()) passOppositeHemi++;
    if (branch.passLargeDeltaEta()) passDeltaEta++;
    if (branch.passLargeDijetMass()) passMjj++;
    if (branch.passAllJetsPtSum()) passPtSum++;
    if (branch.passForwardJet()) passForward++;
    if (branch.passFSRCompatibility()) passFSR++;
    if (branch.passVBFSignalRegion()) passFullSelection++;

    // Print details for first few passing events
    if (branch.passVBFSignalRegion() && passFullSelection <= 5) {
      std::cout << "Event " << totalEvents << " passes VBF selection:" << std::endl;
      std::cout << "  Number of jets: " << branch.countJetsAbovePt(25.0) << std::endl;
      std::cout << "  Leading jet pT: " << branch.topoJetPt[0] << " GeV" << std::endl;
      std::cout << "  Subleading jet pT: " << branch.topoJetPt[1] << " GeV" << std::endl;
      std::cout << "  Dijet mass: " << branch.calcDijetMass() << " GeV" << std::endl;
      std::cout << "  Delta eta_jj: " << std::abs(branch.topoJetEta[0] - branch.topoJetEta[1]) << std::endl;
      std::cout << "  All jets pT sum: " << branch.calcAllJetsPtSum() << " GeV" << std::endl;
      std::cout << std::endl;
    }
  }

  // Print summary
  std::cout << "\n========================================" << std::endl;
  std::cout << "VBF Selection Summary (N=" << totalEvents << " events)" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "Jet multiplicity (2-4, pT>25):  " << passMultiplicity
            << " (" << 100.0*passMultiplicity/totalEvents << "%)" << std::endl;
  std::cout << "Leading jet pT (>60,>40):       " << passLeadPt
            << " (" << 100.0*passLeadPt/totalEvents << "%)" << std::endl;
  std::cout << "Delta phi_jj (<2):              " << passDeltaPhi
            << " (" << 100.0*passDeltaPhi/totalEvents << "%)" << std::endl;
  std::cout << "Opposite hemispheres:           " << passOppositeHemi
            << " (" << 100.0*passOppositeHemi/totalEvents << "%)" << std::endl;
  std::cout << "Delta eta_jj (>3.0):            " << passDeltaEta
            << " (" << 100.0*passDeltaEta/totalEvents << "%)" << std::endl;
  std::cout << "Dijet mass (>600 GeV):          " << passMjj
            << " (" << 100.0*passMjj/totalEvents << "%)" << std::endl;
  std::cout << "All jets pT sum (>140):         " << passPtSum
            << " (" << 100.0*passPtSum/totalEvents << "%)" << std::endl;
  std::cout << "Forward jet (2.38<|eta|<4.0):   " << passForward
            << " (" << 100.0*passForward/totalEvents << "%)" << std::endl;
  std::cout << "FSR compatibility (C_i,m_rel):  " << passFSR
            << " (" << 100.0*passFSR/totalEvents << "%)" << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "FULL VBF SELECTION:             " << passFullSelection
            << " (" << 100.0*passFullSelection/totalEvents << "%)" << std::endl;
  std::cout << "========================================\n" << std::endl;
}

int main(int argc, char** argv) {
  std::string fileNumber = "000001";
  int maxEvents = 1000;

  if (argc > 1) {
    fileNumber = argv[1];
  }
  if (argc > 2) {
    maxEvents = std::atoi(argv[2]);
  }

  test_vbf_selection(fileNumber, maxEvents);
  return 0;
}

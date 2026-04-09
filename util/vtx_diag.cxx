#include <TChain.h>
#include <TTreeReader.h>
#include <iostream>
#include <cmath>
#include <boost/filesystem.hpp>
#include "clustering_includes.h"
#include "clustering_structs.h"
using namespace MyUtl;
int main() {
  TChain chain("ntuple");
  for (const auto& e : boost::filesystem::directory_iterator("../../ntuple-hgtd"))
    chain.Add(e.path().c_str());
  TTreeReader reader(&chain);
  MyUtl::BranchPointerWrapper branch(reader);
  long n=0, n01=0, n02=0, n10=0, n20=0;
  double sum=0;
  while(reader.Next()){
    double dz = std::abs(branch.recoVtxZ[0] - branch.truthVtxZ[0]);
    n++; sum+=dz;
    if(dz < 0.1) n01++;
    if(dz < 0.2) n02++;
    if(dz < 1.0) n10++;
    if(dz < 2.0) n20++;
  }
  std::cout << "Total events: " << n << std::endl;
  std::cout << "|dz| < 0.1mm: " << n01 << "  (" << 100.*n01/n << "%)" << std::endl;
  std::cout << "|dz| < 0.2mm: " << n02 << "  (" << 100.*n02/n << "%)" << std::endl;
  std::cout << "|dz| < 1.0mm: " << n10 << "  (" << 100.*n10/n << "%)" << std::endl;
  std::cout << "|dz| < 2.0mm: " << n20 << "  (" << 100.*n20/n << "%)" << std::endl;
  std::cout << "mean |dz|: " << sum/n << " mm" << std::endl;
  return 0;
}

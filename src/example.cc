// C++ standard library headers
#include <iostream>
#include <vector>

// boost headers
#include <boost/optional.hpp>

// ROOT headers
#include <TFile.h>
#include <TChain.h>
#include <TH1.h>
#include <TLorentzVector.h>

// use <> for library headers, but "" for your own
#include "timed_counter.hh"

using std::cout;
using std::cerr;
using std::endl;

#define TEST(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

int main(int argc, char* argv[]) {
  if (argc<3) { // check program arguments
    cout << "usage: " << argv[0] << " output.root ntuple.root ..." << endl;
    return 1; // returning 1 from main is used to signal error
  }

  // chain input files together
  TChain chain("t3");
  // t3 is the name of the TTree inside the input ROOT files
  // loop over input file names and add them to the chain
  cout << "Input files:" << endl;
  for (int i=2; i<argc; ++i) {
    cout <<"  "<< argv[i] << endl;
    if (!chain.Add(argv[i],0)) return 1; // quit if can't add
  }

  // Define variables
  constexpr size_t Nmax = 4; // maximum number of particles
                             // needs to be larger than max "nparticle"
                             // no penalty if it's a bit larger
  int nparticle, kf[Nmax];
  float px[Nmax], py[Nmax], pz[Nmax], E[Nmax];
  double weight;

  // connect variables to TTree branches
  chain.SetBranchAddress("nparticle", &nparticle); // num particles in event
  chain.SetBranchAddress("kf", kf); // particle MC PID
  chain.SetBranchAddress("E" , E ); // 4-momentum components
  chain.SetBranchAddress("px", px); // no need for &
  chain.SetBranchAddress("py", py); // because float[] decays to float*
  chain.SetBranchAddress("pz", pz);
  chain.SetBranchAddress("weight2", &weight); // event weight
                                              // need to use &

  // open output file
  TFile fout(argv[1],"recreate");
  if (fout.IsZombie()) return 1;
  cout << "Output file: " << fout.GetName() << endl;

  // create histograms
  // name, title, nbins, xmin, xmax
  TH1D *h_H_pT = new TH1D("H_pT","",100,0,1.5e3);
  TH1D *h_Njets_excl = new TH1D("Njets_excl","",Nmax+1,-0.5,0.5+Nmax);
  TH1D *h_Njets_incl = new TH1D("Njets_incl","",Nmax+1,-0.5,0.5+Nmax);
  TH1D *h_jet_pT[Nmax] = {
    new TH1D("jet1_pT","",100,0,1.5e3),
    new TH1D("jet2_pT","",100,0,1.5e3),
    new TH1D("jet3_pT","",100,0,1.5e3),
    new TH1D("jet4_pT","",100,0,1.5e3)
  };

  // particle containers
  boost::optional<TLorentzVector> higgs; // higgs
  std::vector<TLorentzVector> jets; // jets

  // Long64_t is ROOT's typedef for long int
  const Long64_t nentries = chain.GetEntries();
  TEST(nentries)

  using counter = ivanp::timed_counter<Long64_t>;
  for (counter ent(nentries); !!ent; ++ent) {
    chain.GetEntry(ent); // get current TTree entry

    // clear containers
    higgs = boost::none;
    jets.clear();

    // properly assign particles
    for (int i=0; i<nparticle; ++i) {
      if (kf[i] == 25) {
        higgs.emplace(px[i],py[i],pz[i],E[i]);
      } else {
        jets.emplace_back(px[i],py[i],pz[i],E[i]);
      }
    }
    if (!higgs) {
      // print in red
      cerr << "\033[31mNo Higgs in entry " << ent <<"\033[0m"<< endl;
      continue; // move on to next event
    }

    // Fill histograms ==============================================

    const double H_pT = higgs->Pt();
    h_H_pT->Fill(H_pT,weight);

    unsigned Njets = 0; // number of jets that pass cuts

    for (unsigned j=0, n=jets.size(); j<n; ++j) {
      const double jet_pT  = jets[j].Pt();
      const double jet_eta = jets[j].Eta();

      // apply jet cuts
      if (jet_pT  < 30.) continue; // pT cut
      if (jet_eta > 4.4) continue; // eta (pseudo-rapidity) cut

      ++Njets;
      h_jet_pT[j]->Fill(jet_pT,weight);
    }

    h_Njets_excl->Fill(Njets,weight);
    // h_Njets_incl is integral of h_Njets_excl from N to 0
    for (unsigned i=Njets+1; i; ) h_Njets_incl->Fill(--i,weight);

    // ==============================================================
  }

  // fix number of entries for h_Njets_incl
  h_Njets_incl->SetEntries(h_Njets_excl->GetEntries());

  fout.Write(); // write output file
  // Even though the histograms were dynamically allocated,
  // they don't need to be explicitly deleted.
  // The distructor of the TFile they belong to will deallocate them.
}


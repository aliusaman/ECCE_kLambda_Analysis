//____________________________________________________________________________..
//
// This is a template for a Fun4All SubsysReco module with all methods from the
// $OFFLINE_MAIN/include/fun4all/SubsysReco.h baseclass
// You do not have to implement all of them, you can just remove unused methods
// here and in ECCE_DEMP.h.
//
// ECCE_DEMP(const std::string &name = "ECCE_DEMP")
// everything is keyed to ECCE_DEMP, duplicate names do work but it makes
// e.g. finding culprits in logs difficult or getting a pointer to the module
// from the command line
//
// ECCE_kLambda::~ECCE_kLambda()
// this is called when the Fun4AllServer is deleted at the end of running. Be
// mindful what you delete - you do loose ownership of object you put on the node tree
//
// int ECCE_kLambda::Init(PHCompositeNode *topNode)
// This method is called when the module is registered with the Fun4AllServer. You
// can create historgrams here or put objects on the node tree but be aware that
// modules which haven't been registered yet did not put antyhing on the node tree
//
// int ECCE_kLambda::InitRun(PHCompositeNode *topNode)
// This method is called when the first event is read (or generated). At
// this point the run number is known (which is mainly interesting for raw data
// processing). Also all objects are on the node tree in case your module's action
// depends on what else is around. Last chance to put nodes under the DST Node
// We mix events during readback if branches are added after the first event
//
// int ECCE_kLambda::process_event(PHCompositeNode *topNode)
// called for every event. Return codes trigger actions, you find them in
// $OFFLINE_MAIN/include/fun4all/Fun4AllReturnCodes.h
//   everything is good:
//     return Fun4AllReturnCodes::EVENT_OK
//   abort event reconstruction, clear everything and process next event:
//     return Fun4AllReturnCodes::ABORT_EVENT; 
//   proceed but do not save this event in output (needs output manager setting):
//     return Fun4AllReturnCodes::DISCARD_EVENT; 
//   abort processing:
//     return Fun4AllReturnCodes::ABORT_RUN
// all other integers will lead to an error and abort of processing
//
// int ECCE_kLambda::ResetEvent(PHCompositeNode *topNode)
// If you have internal data structures (arrays, stl containers) which needs clearing
// after each event, this is the place to do that. The nodes under the DST node are cleared
// by the framework
//
// int ECCE_kLambda::EndRun(const int runnumber)
// This method is called at the end of a run when an event from a new run is
// encountered. Useful when analyzing multiple runs (raw data). Also called at
// the end of processing (before the End() method)
//
// int ECCE_kLambda::End(PHCompositeNode *topNode)
// This is called at the end of processing. It needs to be called by the macro
// by Fun4AllServer::End(), so do not forget this in your macro
//
// int ECCE_kLambda::Reset(PHCompositeNode *topNode)
// not really used - it is called before the dtor is called
//
// void ECCE_kLambda::Print(const std::string &what) const
// Called from the command line - useful to print information when you need it
//
//____________________________________________________________________________..

#include "ECCE_kLambda.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>

#include <phool/PHCompositeNode.h>

#include <stdio.h>

#include <fun4all/Fun4AllHistoManager.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>

// G4Hits includes
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

// Track includes
#include <trackbase_historic/SvtxTrackMap.h>

// Jet includes
#include <g4eval/JetEvalStack.h>
#include <g4jets/JetMap.h>

// Cluster includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

/// HEPMC truth includes
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

/// Fun4All includes
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/EicEventHeader.h>

#include <TFile.h>
#include <TNtuple.h>
#include <TH2F.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>

#include <algorithm>
#include <cassert>
#include <sstream>
#include <string>
#include <iostream>
#include <stdexcept>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

using namespace std;

ECCE_kLambda::ECCE_kLambda(const std::string &name, const std::string& filename):
 SubsysReco(name)
 , outfilename(filename)
{
  std::cout << "ECCE_kLambda_example::Diff_Tagg_example(const std::string &name) Calling ctor" << std::endl;

  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(m_RandomGenerator, seed);

}

//____________________________________________________________________________..
ECCE_kLambda::~ECCE_kLambda()
{

  gsl_rng_free(m_RandomGenerator);

  std::cout << "ECCE_kLambda::~ECCE_kLambda() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int ECCE_kLambda::Init(PHCompositeNode *topNode)
{
  hm = new Fun4AllHistoManager(Name());
  // create and register your histos (all types) here
  // TH1 *h1 = new TH1F("h1",....)
  // hm->registerHisto(h1);
  outfile = new TFile(outfilename.c_str(), "RECREATE");

  std::cout << "ECCE_kLambda::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  event_itt = 0;

  // create appropriate directories
  
  /*
  // histograms to test events

  h1_trackCharge_Dist = new TH1F("trackCharge_Dist", "Distribution of Charges of Tracks; q [x 1.60 x 10^{-19}]", 10, -2, 2);
  h1_trackPz_Dist = new TH1F("trackPz_Dist", "Distribution of P_{z} of Tracks; p_{z} [GeV/c]", 100, -100, 100);
  h1_nTracks_Dist = new TH1F("nTracks_Dist", "Number of Tracks in Events; # tracks", 10, 0, 10);
  h2_trackPvsTheta_Dist = new TH2F("tracksPvsTheta", "Distribution of P vs #theta for Tracks; p [GeV/c]; #theta [deg]", 100, 0, 100, 100, -180, 180);*/ 

  gDirectory->mkdir("Lambda_Truth");
  gDirectory->cd("Lambda_Truth");
  h1_lTruth_E = new TH1F("lTruth_E", "#Lambda Truth Energy Distribution; E [GeV]", 200, 0, 40);
  h1_lTruth_P = new TH1F("lTruth_P", "#Lambda Truth Momentum Distribution; p [GeV/c]", 200, 0, 40);
  h1_lTruth_Px = new TH1F("lTruth_Px", "#Lambda Truth P_{x} Distribution; p_{x} [GeV/c]", 200, 0, 40);
  h1_lTruth_Py = new TH1F("lTruth_Py", "#Lambda Truth P_{y} Distribution; p_{y} [GeV/c]", 200, 0, 40);
  h1_lTruth_Pz = new TH1F("lTruth_E", "#Lambda Truth P_{z} Distribution; p_{z} [GeV/c]", 200, 0, 40);
  h1_lTruth_Theta = new TH1F("lTruth_Theta", "#Lambda Truth #theta Distribution; #theta [deg]", 200, 0, 3);
  h1_lTruth_Phi = new TH1F("lTruth_Phi", "#Lambda Truth #phi Distribution; #phi [deg]", 200, -180, 180);
  h1_lTruth_PxPy = new TH2F("lTruth_PxPy", "#Lambda Truth P_{x} vs P_{y} Distribution; p_{x} [GeV/c]; p_{y} [GeV/c]", 200, 0, 40, 200, 0, 40);
  h1_lTruth_ThetaP = new TH2F("lTruth_ThetaP", "#Lambda Truth #theta vs P Distribution; #theta [deg]; p[GeV/c]", 200, 0, 3, 200, 0, 40);
  gDirectory->cd("../");

  gDirectory->mkdir("Neutron_Truth");
  gDirectory->cd("Neutron_Truth");
  h1_nTruth_E = new TH1F("nTruth_E", "Neutron Truth Energy Distribution; E [GeV]", 200, 0, 40);
  h1_nTruth_P = new TH1F("nTruth_P", "Neutron Truth Momentum Distribution; p [GeV/c]", 200, 0, 40);
  h1_nTruth_Px = new TH1F("nTruth_Px", "Neutron Truth P_{x} Distribution; p_{x} [GeV/c]", 200, 0, 40);
  h1_nTruth_Py = new TH1F("nTruth_Py", "Neutron Truth P_{y} Distribution; p_{y} [GeV/c]", 200, 0, 40);
  h1_nTruth_Pz = new TH1F("nTruth_E", "Neutron Truth P_{z} Distribution;  p_{z} [GeV/c]", 200, 0, 40);
  h1_nTruth_Theta = new TH1F("nTruth_Theta", "Neutron Truth #theta Distribution; #theta [deg]", 200, 0, 3);
  h1_nTruth_Phi = new TH1F("nTruth_Phi", "Neutron Truth #phi Distribution; #phi [deg]", 200, -180, 180);
  h1_nTruth_PxPy = new TH2F("nTruth_PxPy", "Neutron Truth P_{x} vs P_{y} Distribution; p_{x} [GeV/c]; p_{x} [GeV/c]", 200, 0, 40, 200, 0, 40);
  h1_nTruth_ThetaP = new TH2F("nTruth_ThetaP", "Neutron Truth #theta vs P Distribution; #theta [deg]; p [GeV/c]", 200, 0, 3, 200, 0, 40);
  gDirectory->cd("../");

  gDirectory->mkdir("Electron_Truth");
  gDirectory->cd("Electron_Truth");
  h1_eTruth_E = new TH1F("eTruth_E", "Scattered Electron Truth Energy Distribution; E [GeV]", 200, 0, 8);
  h1_eTruth_P = new TH1F("eTruth_P", "Scattered Electron Truth Momentum Distribution; p [GeV/c]", 200, 0, 8);
  h1_eTruth_Px = new TH1F("eTruth_Px", "Scattered Electron Truth P_{x} Distribution; p_{x} [GeV/c]", 200 , 0, 8);
  h1_eTruth_Py = new TH1F("eTruth_Py", "Scattered Electron Truth P_{y} Distribution; p_{y} [GeV/c]", 200, 0, 8);
  h1_eTruth_Pz = new TH1F("eTruth_Pz", "Scattered Electron Truth P_{z} Distribution; p_{z} [GeV/c]", 200, 0, 8);
  h1_eTruth_Theta = new TH1F("eTruth_Theta", "Scattered Electron Truth #theta Distribution; #theta [deg]", 200, 120, 160);
  h1_eTruth_Phi = new TH1F("eTruth_Phi", "Scattered Electron Truth #phi Distribution; #phi [deg]", 200, -180, 180);
  h1_eTruth_PxPy = new TH2F("eTruth_PxPy", "Scattered Electron Truth P_{x} vs P_{y} Distribution; p_{x} [GeV/c]; p_{y} [GeV/c]", 200, 0, 8, 200, 0, 8);
  h1_eTruth_ThetaP = new TH2F("eTruth_ThetaP", "Scattered Electron Truth #theta vs P Distribution; #theta [deg]; p [GeV/c]", 200, 120, 160, 200, 0, 8);
  gDirectory->cd("../");

  gDirectory->mkdir("Kaon_Truth");
  gDirectory->cd("Kaon_Truth");
  h1_kTruth_E = new TH1F("kTruth_E", "Kaon Truth Energy Distribution; E [GeV]", 200, 0, 25);
  h1_kTruth_P = new TH1F("kTruth_P", "Kaon Truth Momentum Distribution; p [GeV/c]", 200, 0, 25);
  h1_kTruth_Px = new TH1F("kTruth_Px", "Kaon Truth P_{x} Distribution; p_{x} [GeV/c]", 200, 0, 25);
  h1_kTruth_Py = new TH1F("kTruth_Py", "Kaon Truth P_{y} Distribution; p_{y} [GeV/c]", 200, 0, 25);
  h1_kTruth_Pz = new TH1F("kTruth_Pz", "Kaon Truth P_{z} Distribution; p_{z} [GeV/c]", 200, 0, 25);
  h1_kTruth_Theta = new TH1F("kTruth_Theta", "Kaon Truth #theta Distribution; #theta [deg]", 200, 0, 60);
  h1_kTruth_Phi = new TH1F("kTruth_Phi", "Kaon Truth #phi Distribution; #phi [deg]", 200, -180, 180);
  h1_kTruth_PxPy = new TH2F("kTruth_PxPy", "Kaon Truth P_{x} vs P_{y} Distribution; p_{x} [GeV/c]; p_{y} [GeV/c]", 200, 0, 25, 200, 0, 25);
  h1_kTruth_ThetaP = new TH2F("kTruth_ThetaP", "Kaon Truth #theta vs P Distribution; #theta [deg]; p [GeV/c]", 200, 0, 60, 200, 0, 25);
  gDirectory->cd("../");

  /*gDirectory->mkdir("Pion_Truth");
  gDirectory->cd("Pion_Truth");
  h1_pTruth_E = new TH1F("pTruth_E", "Pion Truth Energy Distribution");
  h1_pTruth_P = new TH1F("pTruth_P", "Pion Truth Momentum Distribution");
  h1_pTruth_Px = new TH1F("pTruth_Px", "Pion Truth P_{x} Distribution");
  h1_pTruth_Py = new TH1F("pTruth_Py", "Pion Truth P_{y} Distribution");
  h1_pTruth_Pz = new TH1F("pTruth_E", "Pion Truth P_{z} Distribution");
  h1_pTruth_Theta = new TH1F("pTruth_Theta", "Pion Truth #theta Distribution");
  h1_pTruth_Phi = new TH1F("pTruth_Phi", "Pion Truth #phi Distribution");
  h1_pTruth_PxPy = new TH2F("pTruth_PxPy", "Pion Truth P_{x} vs P_{y} Distribution");
  h1_pTruth_ThetaP = new TH2F("pTruth_ThetaP", "Pion Truth #theta vs P Distribution");
  gDirectory->cd("../");*/

  gDirectory->mkdir("Photon1_Truth");
  gDirectory->cd("Photon1_Truth");
  h1_g1Truth_E = new TH1F("g1Truth_E", "#gamma_{1} Truth Energy Distribution; E [GeV]", 200, -4, 2);
  h1_g1Truth_P = new TH1F("g1Truth_P", "#gamma_{1} Truth Momentum Distribution; p [GeV/c]", 200, 0, 8);
  h1_g1Truth_Px = new TH1F("g1Truth_Px", "#gamma_{1} Truth P_{x} Distribution; p_{x} [GeV/c]", 200, -8, 8);
  h1_g1Truth_Py = new TH1F("g1Truth_Py", "#gamma_{1} Truth P_{y} Distribution; p_{y} [GeV/c]", 200, -8, 8);
  h1_g1Truth_Pz = new TH1F("g1Truth_E", "#gamma_{1} Truth P_{z} Distribution; p_{z} [GeV/c]", 200, -4, 2);
  h1_g1Truth_Theta = new TH1F("g1Truth_Theta", "#gamma_{1} Truth #theta Distribution; #theta [deg]", 200, 90, 120);
  h1_g1Truth_Phi = new TH1F("g1Truth_Phi", "#gamma_{1} Truth #phi Distribution; #phi [deg]", 200, -180, 180);
  h1_g1Truth_PxPy = new TH2F("g1Truth_PxPy", "#gamma_{1} Truth P_{x} vs P_{y} Distribution; p_{x} [GeV/c]; p_{y} [GeV/c]", 200, -8, 8, 200, -8, 8);
  h1_g1Truth_ThetaP = new TH2F("g1Truth_ThetaP", "#gamma_{1} Truth #theta vs P Distribution; #theta [deg]; p [GeV/c]", 200, 90, 120, 200, 0, 8);
  gDirectory->cd("../");

  gDirectory->mkdir("Photon2_Truth");
  gDirectory->cd("Photon2_Truth");
  h1_g2Truth_E = new TH1F("g2Truth_E", "#gamma_{2} Truth Energy Distribution; E [GeV]", 200, -4, 2);
  h1_g2Truth_P = new TH1F("g2Truth_P", "#gamma_{2} Truth Momentum Distribution; p [GeV/c]", 200, 0, 8);
  h1_g2Truth_Px = new TH1F("g2Truth_Px", "#gamma_{2} Truth P_{x} Distribution; p_{x} [GeV/c]", 200, -8, 8);
  h1_g2Truth_Py = new TH1F("g2Truth_Py", "#gamma_{2} Truth P_{y} Distribution; p_{y} [GeV/c]", 200, -8, 8);
  h1_g2Truth_Pz = new TH1F("g2Truth_E", "#gamma_{2} Truth P_{z} Distribution; p_{z} [GeV/c]", 200, -4, 2);
  h1_g2Truth_Theta = new TH1F("g2Truth_Theta", "#gamma_{2} Truth #theta Distribution; #theta [deg]", 200, 90, 120);
  h1_g2Truth_Phi = new TH1F("g2Truth_Phi", "#gamma_{2} Truth #phi Distribution; #phi [deg]", 200, -180, 180);
  h1_g2Truth_PxPy = new TH2F("g2Truth_PxPy", "#gamma_{2} Truth P_{x} vs P_{y} Distribution; p_{x} [GeV/c]; p_{y} [GeV/c]", 200, -8, 8, 200, -8, 8);
  h1_g2Truth_ThetaP = new TH2F("g2Truth_ThetaP", "#gamma_{2} Truth #theta vs P Distribution; #theta [deg]; p [GeV/c]", 200, 90, 120, 200, 0, 8);
  gDirectory->cd("../");

  gDirectory->mkdir("Kinematics");
  gDirectory->cd("Kinematics");
  h1_Q2_Dist = new TH1F("Q2_Dist", "Q^{2} Distribution; Q^{2} [GeV^{2}]", 200, 0, 25);
  h1_W_Dist = new TH1F("W_Dist", "W Distribution; W} [GeV]", 200, 0, 15);
  h1_t_Dist = new TH1F("t_Dist", "-t Distribution; -t [GeV^{2}]", 200, 0, 8);
  h1_xb_Dist = new TH1F("xb_Dist", "x_{b} Distribution; x_{b}", 200, 0, 0.4);
  h1_xi_Dist = new TH1F("xi_Dist", "#xi Distribution; #xi", 200, 0, 0.6);
  gDirectory->cd("../");

  gDirectory->mkdir("Kinematics_Truth");
  gDirectory->cd("Kinematics_Truth");
  h1_Q2Truth_Dist = new TH1F("Q2Truth_Dist", "#frac{#Delta Q^{2}}{Truth Q^{2}} Distribution; (%)", 100, -50, 50);
  h1_WTruth_Dist = new TH1F("WTruth_Dist", "#frac{#Delta W}{Truth W} Distribution; (%)", 100, -50, 50);
  h1_tTruth_Dist = new TH1F("tTruth_Dist", "#frac{#Delta t}{Truth t} Distribution; (%)", 100, -50, 50);
  h1_xbTruth_Dist = new TH1F("xbTruth_Dist", "#frac{#Delta x_{b}}{Truth x_{b}} Distribution; (%)", 100, -50, 50);
  h1_xiTruth_Dist = new TH1F("xiTruth_Dist", "#frac{#Delta #xi}{Truth #xi} Distribution; (%)", 100, -50, 50);
  gDirectory->cd("../");

  gDirectory->mkdir("Kinematics_Coverage");
  gDirectory->cd("Kinematics_Coverage");
  h2_tTruthvQ2_0_25_Dist = new TH2F("tTruthvQ2_0_25_Dist", "#frac{#Delta t}{Truth t} vs Q^{2} Distribution for 0 GeV^{2} <= Q^{2} < 2.5 GeV^{2}; (%); [GeV^{2}]", 100, -50, 50, 200, 0, 25);
  h2_tTruthvQ2_25_375_Dist = new TH2F("tTruthvQ2_25_375_Dist", "#frac{#Delta t}{Truth t} vs Q^{2} Distribution for 2.5 GeV^{2} <= Q^{2} < 3.75 GeV^{2}; (%); [GeV^{2}]", 100, -50, 50, 200, 0, 25);
  h2_tTruthvQ2_375_5_Dist = new TH2F("tTruthvQ2_375_5_Dist", "#frac{#Delta t}{Truth t} vs Q^{2} Distribution for 3.75 GeV^{2} <= Q^{2} < 5 GeV^{2}; (%); [GeV^{2}]", 100, -50, 50, 200, 0, 25);
  h2_tTruthvQ2_5_625_Dist = new TH2F("tTruthvQ2_5_625_Dist", "#frac{#Delta t}{Truth t} vs Q^{2} Distribution for 5 GeV^{2} <= Q^{2} < 6.25 GeV^{2}; (%); [GeV^{2}]", 100, -50, 50, 200, 0, 25);
  h2_tTruthvQ2_625_75_Dist = new TH2F("tTruthvQ2_625_75_Dist", "#frac{#Delta t}{Truth t} vs Q^{2} Distribution for 6.25 GeV^{2} <= Q^{2} < 7.5 GeV^{2}; (%); [GeV^{2}]", 100, -50, 50, 200, 0, 25);
  h2_tTruthvQ2_75_875_Dist = new TH2F("tTruthvQ2_75_875_Dist", "#frac{#Delta t}{Truth t} vs Q^{2} Distribution for 7.5 GeV^{2} <= Q^{2} < 8.75 GeV^{2}; (%); [GeV^{2}]", 100, -50, 50, 200, 0, 25);
  h2_tTruthvQ2_875_10_Dist = new TH2F("tTruthvQ2_875_10_Dist", "#frac{#Delta t}{Truth t} vs Q^{2} Distribution for 8.75 GeV^{2} <= Q^{2} < 10 GeV^{2}; (%); [GeV^{2}]", 100, -50, 50, 200, 0, 25);
  h2_tTruthvQ2_10_1125_Dist = new TH2F("tTruthvQ2_10_1125_Dist", "#frac{#Delta t}{Truth t} vs Q^{2} Distribution for 10 GeV^{2} <= Q^{2} < 11.25 GeV^{2}; (%); [GeV^{2}]", 100, -50, 50, 200, 0, 25);
  h2_tTruthvQ2_g1125_Dist = new TH2F("tTruthvQ2_g1125_Dist", "#frac{#Delta t}{Truth t} vs Q^{2} Distribution for Q^{2} >= 11.25 GeV^{2}; (%); [GeV^{2}]", 100, -50, 50, 200, 0, 25);
  gDirectory->cd("../");

  gDirectory->mkdir("Kinematics_Analysis");
  gDirectory->cd("Kinematics_Analysis");
  h2_t_ep = new TH2F("t_ep", "t vs ScatElec P; t; P_{e'}", 100, 0, 10, 200, 0, 10);
  h2_t_Q2 = new TH2F("t_Q2", "t vs Q^{2}; t; Q^{2}", 100, 0, 10, 200, 0, 50);
  h2_delta_t_t = new TH2F("delta_t_t", "#Delta t vs t; #Delta t (%); t", 200, -100, 100, 100, 0, 1);
  
  for(Int_t A = 0; A < 7; A++){
    h1_t_Q2[A] = new TH1F(Form("t_Q2_%i", (A+1)), Form("t dist, %i < Q^{2} < %i; t", (5 + (A*5)), 10+(A*5)), 100, 0, 10);
    h1_t_alt_Q2[A] = new TH1F(Form("t_alt_Q2_%i", (A+1)), Form("t (Alternative calculation) dist, %i < Q^{2} < %i; t", (5 + (A*5)), 10+(A*5)), 100, 0, 10);
    h2_delta_t_t_Q2[A] = new TH2F(Form("delta_t_t_Q2_%i", (A+1)), Form("#Delta t vs t, %i < Q^{2} < %i; #Delta t (Percent); t", (5 + (A*5)), 10+(A*5)), 200, -100, 100, 100, 0, 1);
  }
  gDirectory->cd("../");

  h2_ZDC_XY = new TH2F("ZDC_XY", "ZDC XY", 200, -50, 50, 200, -50, 50);
  h2_ZDC_XY_Smeared = new TH2F("ZDC_XY_Smeared", "ZDC XY", 200, -50, 50, 200, -50, 50);

  gDirectory->mkdir("Kaon_Information");
  gDirectory->cd("Kaon_Information");
  h1_k_E = new TH1F("k_E", "Kaon Energy Distribution; E [GeV]", 200, 0, 25);
  h1_k_P = new TH1F("k_P", "Kaon Momentum Distribution; p [GeV/c]", 200, 0, 25);
  h1_k_Px = new TH1F("k_Px", "Kaon P_{x} Distribution; p_{x} [GeV/c]", 200, 0, 25);
  h1_k_Py = new TH1F("k_Py", "Kaon P_{y} Distribution; p_{y} [GeV/c]", 200, 0, 25);
  h1_k_Pz = new TH1F("k_E", "Kaon P_{z} Distribution; p_{z} [GeV/c]", 200, 0, 25);
  h1_k_Theta = new TH1F("k_Theta", "Kaon #theta Distribution; #theta [deg]", 200, 0, 60);
  h1_k_Phi = new TH1F("k_Phi", "Kaon #phi Distribution; #phi [deg]", 200, -180, 180);
  h2_k_PxPy = new TH2F("k_PxPy", "Kaon P_{x} vs P_{y} Distribution; p_{x} [GeV/c]; p_{y} [GeV/c]", 200, 0, 25, 200, 0, 25);
  h2_k_ThetaP = new TH2F("k_ThetaP", "Kaon #theta vs P Distribution; #theta [deg]; p [GeV/c]", 200, 0, 60, 200, 0, 25);
  gDirectory->cd("../");

  gDirectory->mkdir("Scattered_Electron_Information");
  gDirectory->cd("Scattered_Electron_Information");
  h1_e_E = new TH1F("e_E", "Electron Energy Distribution; E [GeV]", 200, 0, 25);
  h1_e_P = new TH1F("e_P", "Electron Momentum Distribution; p [GeV/c]", 200, 0, 25);
  h1_e_Px = new TH1F("e_Px", "Electron P_{x} Distribution; p_{x} [GeV/c]", 200, 0, 25);
  h1_e_Py = new TH1F("e_Py", "Electron P_{y} Distribution; p_{y} [GeV/c]", 200, 0, 25);
  h1_e_Pz = new TH1F("e_E", "Electron P_{z} Distribution; p_{z} [GeV/c]", 200, 0, 25);
  h1_e_Theta = new TH1F("e_Theta", "Electron #theta Distribution; #theta [deg]", 200, 0, 60);
  h1_e_Phi = new TH1F("e_Phi", "Electron #phi Distribution; #phi [deg]", 200, -180, 180);
  h2_e_PxPy = new TH2F("e_PxPy", "Electron P_{x} vs P_{y} Distribution; p_{x} [GeV/c]; p_{y} [GeV/c]", 200, 0, 25, 200, 0, 25);
  h2_e_ThetaP = new TH2F("e_ThetaP", "Electron #theta vs P Distribution; #theta [deg]; p [GeV/c]", 200, 0, 60, 200, 0, 25);
  gDirectory->cd("../");

  gDirectory->mkdir("Kaon_Truth_Comparison");
  gDirectory->cd("Kaon_Truth_Comparison");
  h1_kTruthC_E = new TH1F("kTruthC_E", "Kaon #frac{#Delta E}{Truth E} Distribution; (%)", 100, -50, 50);
  h1_kTruthC_P = new TH1F("kTruthC_P", "Kaon #frac{#Delta P}{Truth P} Distribution; (%)", 100, -50, 50);
  h1_kTruthC_Px = new TH1F("kTruthC_Px", "Kaon #frac{#Delta P_{x}}{Truth P_{x}} Distribution; (%)", 100, -50, 50);
  h1_kTruthC_Py = new TH1F("kTruthC_Py", "Kaon #frac{#Delta P_{y}}{Truth P_{y}} Distribution; (%)", 100, -50, 50);
  h1_kTruthC_Pz = new TH1F("kTruthC_Pz", "Kaon #frac{#Delta P_{z}}{Truth P_{z}} Distribution; (%)", 100, -50, 50);
  h1_kTruthC_Theta = new TH1F("kTruthC_Theta", "Kaon #frac{#Delta #theta}{Truth #theta} Distribution; (%)", 100, -50, 50);
  h1_kTruthC_Phi = new TH1F("kTruthC_Phi", "Kaon #frac{#Delta #phi}{Truth #phi} Distribution; (%)", 100, -50, 50);
  h2_kTruthC_PxPy = new TH2F("kTruthC_PxPy", "Kaon #frac{#Delta P_{x}}{Truth P_{x}} vs #frac{#Delta P_{y}}{Truth P_{Y}} Distribution; (%); (%)", 100, -50, 50, 100, -50, 50);
  h2_kTruthC_ThetaP = new TH2F("kTruthC_ThetaP", "Kaon #frac{#Delta #theta}{Truth #theta} vs #frac{#Delta P}{Truth P} Distribution; (%); (%)", 100, -50, 50, 100, -50, 50);
  gDirectory->cd("../");

  gDirectory->mkdir("Electron_Truth_Comparison");
  gDirectory->cd("Electron_Truth_Comparison");
  h1_eTruthC_E = new TH1F("eTruthC_E", "Electron #frac{#Delta E}{Truth E} Distribution; (%)", 100, -50, 50);
  h1_eTruthC_P = new TH1F("eTruthC_P", "Electron #frac{#Delta P}{Truth P} Distribution; (%)", 100, -50, 50);
  h1_eTruthC_Px = new TH1F("eTruthC_Px", "Electron #frac{#Delta P_{x}}{Truth P_{x}} Distribution; (%)", 100, -50, 50);
  h1_eTruthC_Py = new TH1F("eTruthC_Py", "Electron #frac{#Delta P_{y}}{Truth P_{y}} Distribution; (%)", 100, -50, 50);
  h1_eTruthC_Pz = new TH1F("eTruthC_Pz", "Electron #frac{#Delta P_{z}}{Truth P_{z}} Distribution; (%)", 100, -50, 50);
  h1_eTruthC_Theta = new TH1F("eTruthC_Theta", "Electron #frac{#Delta #theta}{Truth #theta} Distribution; (%)", 100, -50, 50);
  h1_eTruthC_Phi = new TH1F("eTruthC_Phi", "Electron #frac{#Delta #phi}{Truth #phi} Distribution; (%)", 100, -50, 50);
  h2_eTruthC_PxPy = new TH2F("eTruthC_PxPy", "Electron #frac{#Delta P_{x}}{Truth P_{x}} vs #frac{#Delta P_{y}}{Truth P_{Y}} Distribution; (%); (%)", 100, -50, 50, 100, -50, 50);
  h2_eTruthC_ThetaP = new TH2F("eTruthC_ThetaP", "Electron #frac{#Delta #theta}{Truth #theta} vs #frac{#Delta P}{Truth P} Distribution; (%); (%)", 100, -50, 50, 100, -50, 50);
  gDirectory->cd("../");

  // Make the histograms which include the weight
  gDirectory->mkdir("Weighted_Distributions");
  gDirectory->cd("Weighted_Distributions");
  gDirectory->mkdir("Kaon_Info");
  gDirectory->cd("Kaon_Info");
  h1_K_px_Weighted = new TH1F("K_px_Weighted", "#K p_{x} Distribution", 200, -20, 20);
  h1_K_py_Weighted = new TH1F("K_py_Weighted", "#K p_{y} Distribution", 200, -20, 20);
  h1_K_pz_Weighted = new TH1F("K_pz_Weighted", "#K p_{z} Distribution", 200, -50, 50); 
  h1_K_p_Weighted = new TH1F("K_p_Weighted", "#K p Distribution", 200, 0, 50);
  h1_K_E_Weighted = new TH1F("K_E_Weighted", "#K E Distribution", 200, 0, 50);
  h1_K_Theta_Weighted = new TH1F("K_Theta_Weighted", "#K #theta Distribution; #theta [deg]", 200, 0, 50);
  h1_K_Phi_Weighted = new TH1F("K_Phi_Weighted", "#K #phi Distribution; #phi [deg]", 360, -180, 180);
  h2_KTrack_ThetaPhi_Weighted = new TH2F("KTrack_ThetaPhi_Weighted", "#K Track #theta vs #phi; #theta [deg]; #phi [deg]", 120, 0, 60, 720, -180, 180);
  h2_KTrack_pTheta_Weighted = new TH2F("KTrack_pTheta_Weighted", "#K Track #theta vs P; #theta [deg]; P [GeV/c]", 120, 0, 60, 500, 0, 50);
  h2_KTrack_ThetaPhi_Smeared_Weighted = new TH2F("KTrack_ThetaPhi_Smeared_Weighted", "#K Track #theta vs #phi; #theta [deg]; #phi [deg]", 120, 0, 60, 720, -180, 180);
  h2_KTrack_pTheta_Smeared_Weighted = new TH2F("KTrack_pTheta_Smeared_Weighted", "#K Track #theta vs P; #theta [deg]; P [GeV/c]", 120, 0, 60, 500, 0, 50);
  gDirectory->cd("../");

  gDirectory->mkdir("Kaon_Truth_Info");
  gDirectory->cd("Kaon_Truth_Info");
  h1_KTruth_p_Weighted = new TH1F("KTruth_p_Weighted", "#K #frac{#Delta p}{Truth p} Distribution (%); %", 100, -50, 50);
  h1_KTruth_px_Weighted = new TH1F("KTruth_px_Weighted", "#K #frac{#Delta px}{Truth px} Distribution (%); %", 100, -50, 50);
  h1_KTruth_py_Weighted = new TH1F("KTruth_py_Weighted", "#K #frac{#Delta py}{Truth py} Distribution (%); %", 100, -50, 50);
  h1_KTruth_pz_Weighted = new TH1F("KTruth_pz_Weighted", "#K #frac{#Delta pz}{Truth pz} Distribution (%); %", 100, -50, 50);
  h1_KTruth_E_Weighted = new TH1F("KTruth_E_Weighted", "#K #frac{#Delta E}{Truth E} Distribution (%); %", 100, -50, 50);
  h1_KTruth_p_Smeared_Weighted = new TH1F("KTruth_p_Smeared_Weighted", "#K #frac{#Delta p}{Truth p} Distribution (%); %", 100, -50, 50);
  h1_KTruth_px_Smeared_Weighted = new TH1F("KTruth_px_Smeared_Weighted", "#K #frac{#Delta px}{Truth px} Distribution (%); %", 100, -50, 50);
  h1_KTruth_py_Smeared_Weighted = new TH1F("KTruth_py_Smeared_Weighted", "#K #frac{#Delta py}{Truth py} Distribution (%); %", 100, -50, 50);
  h1_KTruth_pz_Smeared_Weighted = new TH1F("KTruth_pz_Smeared_Weighted", "#K #frac{#Delta pz}{Truth pz} Distribution (%); %", 100, -50, 50);
  h1_KTruth_E_Smeared_Weighted = new TH1F("KTruth_E_Smeared_Weighted", "#K #frac{#Delta E}{Truth E} Distribution (%); %", 100, -50, 50);
  h2_KTruth_pxpy_Weighted = new TH2F("KTruth_pxpy_Weighted", "#K #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);
  h2_KTruth_pxpz_Weighted = new TH2F("KTruth_pxpz_Weighted", "#K #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);
  h2_KTruth_pypz_Weighted = new TH2F("KTruth_pypz_Weighted", "#K #frac{#Delta p_{y}}{Truth p_{y}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);
  h2_KTruth_pxpy_Smeared_Weighted = new TH2F("KTruth_pxpy_Smeared_Weighted", "#K #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);
  gDirectory->cd("../");
  
  gDirectory->mkdir("Electron_Info");
  gDirectory->cd("Scattered_Electron_Info");
  h1_e_px_Weighted = new TH1F("e_px_Weighted", "e' p_{x} Distribution", 200, -10, 10);
  h1_e_py_Weighted = new TH1F("e_py_Weighted", "e' p_{y} Distribution", 200, -10, 10);
  h1_e_pz_Weighted = new TH1F("e_pz_Weighted", "e' p_{z} Distribution", 200, -10, 0); 
  h1_e_p_Weighted = new TH1F("e_p_Weighted", "e' p Distribution", 200, 0, 10);
  h1_e_E_Weighted = new TH1F("e_E_Weighted", "e' E Distribution", 200, 0, 10);
  h1_e_Theta_Weighted = new TH1F("e_Theta_Weighted", "e' #theta Distribution; #theta [deg]", 200, 110, 160);
  h1_e_Phi_Weighted = new TH1F("e_Phi_Weighted", "e' #phi Distribution; #phi [deg]", 360, -180, 180);
  h2_eTrack_ThetaPhi_Weighted = new TH2F("eTrack_ThetaPhi_Weighted", "e' Track #theta vs #phi; #theta [deg]; #phi [deg]", 140, 110, 180, 720, -180, 180);
  h2_eTrack_pTheta_Weighted = new TH2F("eTrack_pTheta_Weighted", "e' Track #theta vs P; #theta [deg]; P [GeV/c]", 140, 110, 180, 100, 0, 10);
  h2_eTrack_ThetaPhi_Smeared_Weighted = new TH2F("eTrack_ThetaPhi_Smeared_Weighted", "e' Track #theta vs #phi; #theta [deg]; #phi [deg]", 140, 110, 180, 720, -180, 180);
  h2_eTrack_pTheta_Smeared_Weighted = new TH2F("eTrack_pTheta_Smeared_Weighted", "e' Track #theta vs P; #theta [deg]; P [GeV/c]", 140, 110, 180, 100, 0, 10);
  gDirectory->cd("../");

  gDirectory->mkdir("Electron_Truth_Info");
  gDirectory->cd("Scattered_Electron_Truth_Info");
  h1_eTruth_p_Weighted = new TH1F("eTruth_p_Weighted", "e' #frac{#Delta p}{Truth p} Distribution (%) ; %", 100, -50, 50);
  h1_eTruth_px_Weighted = new TH1F("eTruth_px_Weighted", "#e' #frac{#Delta px}{Truth px} Distribution (%); %", 100, -50, 50);
  h1_eTruth_py_Weighted = new TH1F("eTruth_py_Weighted", "#e' #frac{#Delta py}{Truth py} Distribution (%); %", 100, -50, 50);
  h1_eTruth_pz_Weighted = new TH1F("eTruth_pz_Weighted", "e' #frac{#Delta pz}{Truth pz} Distribution (%); %", 100, -50, 50);
  h1_eTruth_E_Weighted = new TH1F("eTruth_E_Weighted", "e' #frac{#Delta E}{Truth E} Distribution (%) ; %", 100, -50, 50);
  h1_eTruth_p_Smeared_Weighted = new TH1F("eTruth_p_Smeared_Weighted", "e' #frac{#Delta p}{Truth p} Distribution (%) ; %", 100, -50, 50);
  h1_eTruth_px_Smeared_Weighted = new TH1F("eTruth_px_Smeared_Weighted", "#e' #frac{#Delta px}{Truth px} Distribution (%); %", 100, -50, 50);
  h1_eTruth_py_Smeared_Weighted = new TH1F("eTruth_py_Smeared_Weighted", "#e' #frac{#Delta py}{Truth py} Distribution (%); %", 100, -50, 50);
  h1_eTruth_pz_Smeared_Weighted = new TH1F("eTruth_pz_Smeared_Weighted", "e' #frac{#Delta pz}{Truth pz} Distribution (%); %", 100, -50, 50);
  h1_eTruth_E_Smeared_Weighted = new TH1F("eTruth_E_Smeared_Weighted", "e' #frac{#Delta E}{Truth E} Distribution (%) ; %", 100, -50, 50);
  h2_eTruth_pxpy_Weighted = new TH2F("eTruth_pxpy_Weighted", "e' #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);  
  h2_eTruth_pxpz_Weighted = new TH2F("eTruth_pxpz_Weighted", "e' #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);  
  h2_eTruth_pypz_Weighted = new TH2F("eTruth_pypz_Weighted", "e' #frac{#Delta p_{y}}{Truth p_{y}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);  
  h2_eTruth_pxpy_Smeared_Weighted = new TH2F("eTruth_pxpy_Smeared_Weighted", "e' #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);
  gDirectory->cd("../");

  gDirectory->mkdir("Lambda_Info");
  gDirectory->cd("Lambda_Info");
  h1_L_px_Weighted = new TH1F("L_px_Weighted", "L p_{x} Distribution", 320, -4, 4);
  h1_L_py_Weighted = new TH1F("L_py_Weighted", "L p_{y} Distribution", 200, -2.5, 2.5);
  h1_L_pz_Weighted = new TH1F("L_pz_Weighted", "L p_{z} Distribution", 240, 0, 120); 
  h1_L_p_Weighted = new TH1F("L_p_Weighted", "L p Distribution", 240, 0, 120);
  h1_L_E_Weighted = new TH1F("L_E_Weighted", "L E Distribution", 240, 0, 120);
  h1_L_Theta_Weighted = new TH1F("L_Theta_Weighted", "L #theta Distribution; #theta [deg]", 500, 0, 5);
  h1_L_Phi_Weighted = new TH1F("L_Phi_Weighted", "L #phi Distribution; #phi [deg]", 360, -180, 180);
  h2_LTrack_ThetaPhi_Weighted = new TH2F("LTrack_ThetaPhi_Weighted", "L Track #theta vs #phi; #theta [deg]; #phi [deg]", 500, 0, 5, 360, -180, 180);
  h2_LTrack_pTheta_Weighted = new TH2F("LTrack_pTheta_Weighted", "L Track #theta vs P; #theta [deg]; P [GeV/c]", 100, 0, 1, 1000, 0, 100);
  h2_LTrack_ThetaPhi_Smeared_Weighted = new TH2F("LTrack_ThetaPhi_Smeared_Weighted", "L Track #theta vs #phi; #theta [deg]; #phi [deg]", 100, 0, 1, 100, -50, 50);
  h2_LTrack_pTheta_Smeared_Weighted = new TH2F("LTrack_pTheta_Smeared_Weighted", "L Track #theta vs P; #theta [deg]; P [GeV/c]", 100, 0, 1, 1000, 0, 100);
  gDirectory->cd("../");

  gDirectory->mkdir("Lambda_Truth_Info");
  gDirectory->cd("Lambda_Truth_Info");
  h1_LTruth_p_Weighted = new TH1F("LTruth_p_Weighted", "L #frac{#Delta p}{Truth p} Distribution (%) ; %", 100, -50, 50);
  h1_LTruth_px_Weighted = new TH1F("LTruth_px_Weighted", "#L #frac{#Delta px}{Truth px} Distribution (%); %", 100, -50, 50);
  h1_LTruth_py_Weighted = new TH1F("LTruth_py_Weighted", "#L #frac{#Delta py}{Truth py} Distribution (%); %", 100, -50, 50);
  h1_LTruth_pz_Weighted = new TH1F("LTruth_pz_Weighted", "L #frac{#Delta pz}{Truth pz} Distribution (%); %", 100, -50, 50);
  h1_LTruth_E_Weighted = new TH1F("LTruth_E_Weighted", "L #frac{#Delta E}{Truth E} Distribution (%) ; %", 100, -50, 50);
  h1_LTruth_p_Smeared_Weighted = new TH1F("LTruth_p_Smeared_Weighted", "L #frac{#Delta p}{Truth p} Distribution (%) ; %", 100, -50, 50);
  h1_LTruth_px_Smeared_Weighted = new TH1F("LTruth_px_Smeared_Weighted", "#L #frac{#Delta px}{Truth px} Distribution (%); %", 100, -50, 50);
  h1_LTruth_py_Smeared_Weighted = new TH1F("LTruth_py_Smeared_Weighted", "#L #frac{#Delta py}{Truth py} Distribution (%); %", 100, -50, 50);
  h1_LTruth_pz_Smeared_Weighted = new TH1F("LTruth_pz_Smeared_Weighted", "L #frac{#Delta pz}{Truth pz} Distribution (%); %", 100, -50, 50);
  h1_LTruth_E_Smeared_Weighted = new TH1F("LTruth_E_Smeared_Weighted", "L #frac{#Delta E}{Truth E} Distribution (%) ; %", 100, -50, 50);
  // SJDK 04/08/21 - Neutron Px distributions currently don't work, a correction is applied for the ZDC x position which is used to determine the neutron 4 vector (from angles/E). This is NOT applied to the truth track. Need to figure out how to correctly adjust the truth track too
  h2_LTruth_pxpy_Weighted = new TH2F("LTruth_pxpy_Weighted", "L #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);
  h2_LTruth_pxpz_Weighted = new TH2F("LTruth_pxpz_Weighted", "L #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);
  h2_LTruth_pypz_Weighted = new TH2F("LTruth_pypz_Weighted", "L #frac{#Delta p_{y}}{Truth p_{y}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);
  h2_LTruth_pxpy_Smeared_Weighted = new TH2F("LTruth_pxpy_Smeared_Weighted", "L #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);
  gDirectory->cd("../");

  gDirectory->mkdir("Neutron_Info");
  gDirectory->cd("Neutron_Info");
  h1_n_px_Weighted = new TH1F("n_px_Weighted", "n p_{x} Distribution", 320, -4, 4);
  h1_n_py_Weighted = new TH1F("n_py_Weighted", "n p_{y} Distribution", 200, -2.5, 2.5);
  h1_n_pz_Weighted = new TH1F("n_pz_Weighted", "n p_{z} Distribution", 240, 0, 120); 
  h1_n_p_Weighted = new TH1F("n_p_Weighted", "n p Distribution", 240, 0, 120);
  h1_n_E_Weighted = new TH1F("n_E_Weighted", "n E Distribution", 240, 0, 120);
  h1_n_Theta_Weighted = new TH1F("n_Theta_Weighted", "n #theta Distribution; #theta [deg]", 500, 0, 5);
  h1_n_Phi_Weighted = new TH1F("n_Phi_Weighted", "n #phi Distribution; #phi [deg]", 360, -180, 180);
  h2_nTrack_ThetaPhi_Weighted = new TH2F("nTrack_ThetaPhi_Weighted", "n Track #theta vs #phi; #theta [deg]; #phi [deg]", 500, 0, 5, 360, -180, 180);
  h2_nTrack_pTheta_Weighted = new TH2F("nTrack_pTheta_Weighted", "n Track #theta vs P; #theta [deg]; P [GeV/c]", 100, 0, 1, 1000, 0, 100);
  h2_nTrack_ThetaPhi_Smeared_Weighted = new TH2F("nTrack_ThetaPhi_Smeared_Weighted", "n Track #theta vs #phi; #theta [deg]; #phi [deg]", 100, 0, 1, 100, -50, 50);
  h2_nTrack_pTheta_Smeared_Weighted = new TH2F("nTrack_pTheta_Smeared_Weighted", "n Track #theta vs P; #theta [deg]; P [GeV/c]", 100, 0, 1, 1000, 0, 100);
  gDirectory->cd("../");

  gDirectory->mkdir("Neutron_Truth_Info");
  gDirectory->cd("Neutron_Truth_Info");
  h1_nTruth_p_Weighted = new TH1F("nTruth_p_Weighted", "n #frac{#Delta p}{Truth p} Distribution (%) ; %", 100, -50, 50);
  h1_nTruth_px_Weighted = new TH1F("nTruth_px_Weighted", "#n #frac{#Delta px}{Truth px} Distribution (%); %", 100, -50, 50);
  h1_nTruth_py_Weighted = new TH1F("nTruth_py_Weighted", "#n #frac{#Delta py}{Truth py} Distribution (%); %", 100, -50, 50);
  h1_nTruth_pz_Weighted = new TH1F("nTruth_pz_Weighted", "n #frac{#Delta pz}{Truth pz} Distribution (%); %", 100, -50, 50);
  h1_nTruth_E_Weighted = new TH1F("nTruth_E_Weighted", "n #frac{#Delta E}{Truth E} Distribution (%) ; %", 100, -50, 50);
  h1_nTruth_p_Smeared_Weighted = new TH1F("nTruth_p_Smeared_Weighted", "n #frac{#Delta p}{Truth p} Distribution (%) ; %", 100, -50, 50);
  h1_nTruth_px_Smeared_Weighted = new TH1F("nTruth_px_Smeared_Weighted", "#n #frac{#Delta px}{Truth px} Distribution (%); %", 100, -50, 50);
  h1_nTruth_py_Smeared_Weighted = new TH1F("nTruth_py_Smeared_Weighted", "#n #frac{#Delta py}{Truth py} Distribution (%); %", 100, -50, 50);
  h1_nTruth_pz_Smeared_Weighted = new TH1F("nTruth_pz_Smeared_Weighted", "n #frac{#Delta pz}{Truth pz} Distribution (%); %", 100, -50, 50);
  h1_nTruth_E_Smeared_Weighted = new TH1F("nTruth_E_Smeared_Weighted", "n #frac{#Delta E}{Truth E} Distribution (%) ; %", 100, -50, 50);
  // SJDK 04/08/21 - Neutron Px distributions currently don't work, a correction is applied for the ZDC x position which is used to determine the neutron 4 vector (from angles/E). This is NOT applied to the truth track. Need to figure out how to correctly adjust the truth track too
  h2_nTruth_pxpy_Weighted = new TH2F("nTruth_pxpy_Weighted", "n #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);
  h2_nTruth_pxpz_Weighted = new TH2F("nTruth_pxpz_Weighted", "n #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);
  h2_nTruth_pypz_Weighted = new TH2F("nTruth_pypz_Weighted", "n #frac{#Delta p_{y}}{Truth p_{y}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);
  h2_nTruth_pxpy_Smeared_Weighted = new TH2F("nTruth_pxpy_Smeared_Weighted", "n #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);
  gDirectory->cd("../");

  gDirectory->mkdir("Photon1_Info");
  gDirectory->cd("Photon1_Info");
  h1_g1_px_Weighted = new TH1F("g1_px_Weighted", "g1 p_{x} Distribution", 320, -4, 4);
  h1_g1_py_Weighted = new TH1F("g1_py_Weighted", "g1 p_{y} Distribution", 200, -2.5, 2.5);
  h1_g1_pz_Weighted = new TH1F("g1_pz_Weighted", "g1 p_{z} Distribution", 240, 0, 120); 
  h1_g1_p_Weighted = new TH1F("g1_p_Weighted", "g1 p Distribution", 240, 0, 120);
  h1_g1_E_Weighted = new TH1F("g1_E_Weighted", "g1 E Distribution", 240, 0, 120);
  h1_g1_Theta_Weighted = new TH1F("g1_Theta_Weighted", "g1 #theta Distribution; #theta [deg]", 500, 0, 5);
  h1_g1_Phi_Weighted = new TH1F("g1_Phi_Weighted", "g1 #phi Distribution; #phi [deg]", 360, -180, 180);
  h2_g1Track_ThetaPhi_Weighted = new TH2F("g1Track_ThetaPhi_Weighted", "g1 Track #theta vs #phi; #theta [deg]; #phi [deg]", 500, 0, 5, 360, -180, 180);
  h2_g1Track_pTheta_Weighted = new TH2F("g1Track_pTheta_Weighted", "g1 Track #theta vs P; #theta [deg]; P [GeV/c]", 100, 0, 1, 1000, 0, 100);
  h2_g1Track_ThetaPhi_Smeared_Weighted = new TH2F("g1Track_ThetaPhi_Smeared_Weighted", "g1 Track #theta vs #phi; #theta [deg]; #phi [deg]", 100, 0, 1, 100, -50, 50);
  h2_g1Track_pTheta_Smeared_Weighted = new TH2F("g1Track_pTheta_Smeared_Weighted", "g1 Track #theta vs P; #theta [deg]; P [GeV/c]", 100, 0, 1, 1000, 0, 100);
  gDirectory->cd("../");

  gDirectory->mkdir("Photon2_Truth_Info");
  gDirectory->cd("Photon2_Truth_Info");
  h1_g2Truth_p_Weighted = new TH1F("g2Truth_p_Weighted", "g2 #frac{#Delta p}{Truth p} Distribution (%) ; %", 100, -50, 50);
  h1_g2Truth_px_Weighted = new TH1F("g2Truth_px_Weighted", "#g2 #frac{#Delta px}{Truth px} Distribution (%); %", 100, -50, 50);
  h1_g2Truth_py_Weighted = new TH1F("g2Truth_py_Weighted", "#g2 #frac{#Delta py}{Truth py} Distribution (%); %", 100, -50, 50);
  h1_g2Truth_pz_Weighted = new TH1F("g2Truth_pz_Weighted", "g2 #frac{#Delta pz}{Truth pz} Distribution (%); %", 100, -50, 50);
  h1_g2Truth_E_Weighted = new TH1F("g2Truth_E_Weighted", "g2 #frac{#Delta E}{Truth E} Distribution (%) ; %", 100, -50, 50);
  h1_g2Truth_p_Smeared_Weighted = new TH1F("g2Truth_p_Smeared_Weighted", "g2 #frac{#Delta p}{Truth p} Distribution (%) ; %", 100, -50, 50);
  h1_g2Truth_px_Smeared_Weighted = new TH1F("g2Truth_px_Smeared_Weighted", "#g2 #frac{#Delta px}{Truth px} Distribution (%); %", 100, -50, 50);
  h1_g2Truth_py_Smeared_Weighted = new TH1F("g2Truth_py_Smeared_Weighted", "#g2 #frac{#Delta py}{Truth py} Distribution (%); %", 100, -50, 50);
  h1_g2Truth_pz_Smeared_Weighted = new TH1F("g2Truth_pz_Smeared_Weighted", "g2 #frac{#Delta pz}{Truth pz} Distribution (%); %", 100, -50, 50);
  h1_g2Truth_E_Smeared_Weighted = new TH1F("g2Truth_E_Smeared_Weighted", "g2 #frac{#Delta E}{Truth E} Distribution (%) ; %", 100, -50, 50);
  // SJDK 04/08/21 - Neutron Px distributions currently don't work, a correction is applied for the ZDC x position which is used to determine the neutron 4 vector (from angles/E). This is NOT applied to the truth track. Need to figure out how to correctly adjust the truth track too
  h2_g2Truth_pxpy_Weighted = new TH2F("g2Truth_pxpy_Weighted", "g2 #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);
  h2_g2Truth_pxpz_Weighted = new TH2F("g2Truth_pxpz_Weighted", "g2 #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);
  h2_g2Truth_pypz_Weighted = new TH2F("g2Truth_pypz_Weighted", "g2 #frac{#Delta p_{y}}{Truth p_{y}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);
  h2_g2Truth_pxpy_Smeared_Weighted = new TH2F("g2Truth_pxpy_Smeared_Weighted", "g2 #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);
  gDirectory->cd("../");

  gDirectory->mkdir("Kinematics_Info");
  gDirectory->cd("Kinematics_Info");
  h1_Q2_Dist_Weighted = new TH1F("Q2_Dist_Weighted", "Q^{2} Distribution", 200, 0, 50);
  h1_W_Dist_Weighted = new TH1F("W_Dist_Weighted", "W Distribution", 500, 0, 50);
  h1_t_Dist_Weighted = new TH1F("t_Dist_Weighted", "t Distribution", 100, 0, 10);
  h1_t_alt_Dist_Weighted = new TH1F("t_alt_Dist_Weighted", "t (Alternative calculation) Distribution", 100, 0, 10);
  h1_t_comp_Weighted = new TH1F("t_comp_Dist_Weighted", "#frac{#Delta t}{t} Distribution; #frac{t_{alt}-t}{t} (%)", 200, -100, 100);
  h1_xb_Dist_Weighted = new TH1F("xb_Dist_Weighted", "x_{b} Distribution", 100, 0, 1);
  h1_xi_Dist_Weighted = new TH1F("xi_Dist_Weighted", "#xi Distribution", 100, 0, 1);
  gDirectory->cd("../");

  gDirectory->mkdir("Kinematics_Truth_Info");
  gDirectory->cd("Kinematics_Truth_Info");
  h1_Q2Truth_Dist_Weighted = new TH1F("Q2Truth_Dist_Weighted", "Q^{2} Truth Distribution", 200, 0, 50);
  h1_WTruth_Dist_Weighted = new TH1F("WTruth_Dist_Weighted", "W Truth Distribution", 500, 0, 50);
  h1_tTruth_Dist_Weighted = new TH1F("tTruth_Dist_Weighted", "t Truth Distribution", 100, 0, 10);
  h1_xbTruth_Dist_Weighted = new TH1F("xbTruth_Dist_Weighted", "x_{b} Truth Distribution", 100, 0, 1);
  h1_xiTruth_Dist_Weighted = new TH1F("xiTruth_Dist_Weighted", "#xi Truth Distribution", 100, 0, 1);
  gDirectory->cd("../");

  gDirectory->mkdir("Kinematics_Analysis");
  gDirectory->cd("Kinematics_Analysis");
  h2_t_ep_Weighted = new TH2F("t_ep_Weighted", "t vs ScatElec P; t; P_{e'}", 100, 0, 10, 200, 0, 10);
  h2_t_Q2_Weighted = new TH2F("t_Q2_Weighted", "t vs Q^{2}; t; Q^{2}", 100, 0, 10, 200, 0, 50);
  h2_delta_t_t_Weighted = new TH2F("delta_t_t_Weighted", "#Delta t vs t; #Delta t (%); t", 200, -100, 100, 100, 0, 1);
  
  for(Int_t A = 0; A < 7; A++){
    h1_t_Q2_Weighted[A] = new TH1F(Form("t_Q2_%i_Weighted", (A+1)), Form("t dist, %i < Q^{2} < %i; t", (5 + (A*5)), 10+(A*5)), 100, 0, 10);
    h1_t_alt_Q2_Weighted[A] = new TH1F(Form("t_alt_Q2_%i_Weighted", (A+1)), Form("t (Alternative calculation) dist, %i < Q^{2} < %i; t", (5 + (A*5)), 10+(A*5)), 100, 0, 10);
    h2_delta_t_t_Q2_Weighted[A] = new TH2F(Form("delta_t_t_Q2_%i_Weighted", (A+1)), Form("#Delta t vs t, %i < Q^{2} < %i; #Delta t (Percent); t", (5 + (A*5)), 10+(A*5)), 200, -100, 100, 100, 0, 1);
  }
  gDirectory->cd("../");

  h2_ZDC_XY_Weighted = new TH2F("ZDC_XY_Weighted", "ZDC XY", 200, -50, 50, 200, -50, 50);
  h2_ZDC_XY_Smeared_Weighted = new TH2F("ZDC_XY_Smeared_Weighted", "ZDC XY", 200, -50, 50, 200, -50, 50);

  // Define beam 4 vectors
  e_beam_energy = 5;
  e_beam_pmag = sqrt(pow(e_beam_energy,2)-pow(mElec,2));
  ion_beam_energy = 41;
  ion_beam_pmag = sqrt((pow(ion_beam_energy,2)-pow(mProt,2)));
  crossing_angle = 0.05; 
  Double_t Pi = TMath::ACos(-1);
  eBeam4Vect.SetPxPyPzE(0,0,-1*e_beam_pmag,e_beam_energy);
  pBeam4Vect.SetPxPyPzE(-ion_beam_pmag*TMath::Sin(crossing_angle),ion_beam_pmag*TMath::Sin(crossing_angle)*TMath::Sin(Pi),ion_beam_pmag*TMath::Cos(crossing_angle),ion_beam_energy);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_kLambda::InitRun(PHCompositeNode *topNode)
{
  std::cout << "ECCE_kLambda::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_kLambda::process_event(PHCompositeNode *topNode)
{
  EEMC_hit = 0;
  //  Double_t Pi = TMath::ACos(-1);
  //  eBeam4Vect.SetPxPyPzE(0,0,-5,5);
  //  pBeam4Vect.SetPxPyPzE(-41*TMath::Sin(0.05),41*TMath::Sin(0.05)*TMath::Sin(Pi),41*TMath::Cos(0.05),41); // 5 on 41

  event_itt++; 
 
  if(event_itt%100 == 0)
     std::cout << "Event Processing Counter: " << event_itt << endl;

  // made considering there is no ZDC


  /*
  // event testing code
  // Get track map

  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  if (!trackmap)
    {
      trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
      if (!trackmap)
    	{
    	  cout << "ECCE_kLambda::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
    	  exit(-1);
    	}
    }

  Int_t numberTracks = 0;

  // Loop over our tracks, assign info
  for (SvtxTrackMap::Iter iter = trackmap->begin(); iter != trackmap->end();++iter) {
    SvtxTrack* track = iter->second;
    h1_trackCharge_Dist->Fill(track->get_charge());
    h1_trackPz_Dist->Fill(track->get_pz());
    trackVect.SetXYZ(track->get_px(), track->get_py(), track->get_pz());
    track4Vect.SetPxPyPzE(track->get_px(), track->get_py(), track->get_pz(), sqrt(pow(trackVect.Mag(), 2)+pow(mKaonP,2)));
    h2_trackPvsTheta_Dist->Fill(track4Vect.P(),track4Vect.Theta()*TMath::RadToDeg());
    numberTracks++;
    }

  h1_nTracks_Dist->Fill(numberTracks);
  */

  if (Check_hits(topNode) == true && Check_eKaon(topNode) == true){
    // Get track map
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
    // Get MC truth info
    PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    // Get ZDC hits
    //    PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_ZDC");
    PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_ZDCsurrogate");
    // Get the primary particle range
    PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
    //Get the secondary particle range
    PHG4TruthInfoContainer::Range range_s = truthinfo->GetSecondaryParticleRange();

    if (!trackmap) {
      trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
      if (!trackmap){
    	cout << "ECCE_kLambda::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
    	exit(-1);
      }
    }

    if (!truthinfo) {
      cout << PHWHERE << "PHG4TruthInfoContainer node is missing, can't collect G4 truth particles" << endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }

    /* if (hits) {
      std::cout << "zdc" << endl;
    } */

    /// Loop over the G4 truth (stable) particles
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
      {
	/// Get this truth particle
	const PHG4Particle *truth = iter->second;
	if ( truth->get_pid() == 11){ // PDG 11 -> Scattered electron
	  e4VectTruth.SetPxPyPzE(truth->get_px(), truth->get_py(), truth->get_pz(), truth->get_e());
	}
	else if (truth->get_pid() == 321){ // PDG 321 -> Kaon^+
	  kaon4VectTruth.SetPxPyPzE(truth->get_px(), truth->get_py(), truth->get_pz(), truth->get_e());
	}
	else if (truth->get_pid() == 3122){ // PDG 3122 -> Lambda
	  l4VectTruth.SetPxPyPzE(truth->get_px(), truth->get_py(), truth->get_pz(), truth->get_e());
	}
      }

    /// Determine if there are daughter tracks

    Double_t ZDCEnergy = 0;
    Bool_t ZDCHits = kFALSE;
    Int_t nHits = 0;
    Int_t Gamma = kFALSE;

    for (PHG4TruthInfoContainer::ConstIterator iter = range.first;iter != range.second;++iter)
    {
	/// Get this truth particle
	const PHG4Particle *truth = iter->second;

        if (truth->get_pid() == 3122){ // PDG 3122 -> Lambda
          for (PHG4TruthInfoContainer::ConstIterator iter_s = range_s.first; iter_s !=range_s.second; ++iter_s) {
	    const PHG4Particle *truth_s = iter_s->second;

             if ( truth->get_pid() == 2112){ // PDG 2112 -> Neutron
  	       n4VectTruth.SetPxPyPzE(truth_s->get_px(), truth_s->get_py(), truth_s->get_pz(), truth_s->get_e());
	       nHits++;
	       ZDCEnergy = ZDCEnergy + n4VectTruth.E();
             }

             else if ( truth->get_pid() == 22 && Gamma == kTRUE){ // PDG 22 -> Photon (number 2)
	       g24VectTruth.SetPxPyPzE(truth_s->get_px(), truth_s->get_py(), truth_s->get_pz(), truth_s->get_e());
	       nHits++;
	       ZDCEnergy = ZDCEnergy + g24VectTruth.E();
	     }

             else if ( truth->get_pid() == 22){ // PDG 22 -> Photon (number 1)
               g14VectTruth.SetPxPyPzE(truth_s->get_px(), truth_s->get_py(), truth_s->get_pz(), truth_s->get_e());
	       nHits++;
	       ZDCEnergy = ZDCEnergy + g14VectTruth.E();
	       Gamma = kTRUE;
	     }
	  }

	  Bool_t Energy = kFALSE;

	  if (nHits == 3 && ZDCEnergy >= 10) {
	    Energy = kTRUE;
	  }

	  else if (l4VectTruth.E() >= 10) {
	    Energy = kTRUE;
          }
	      
          Double_t diffPhi = TMath::Abs((g14VectTruth.Phi()*TMath::RadToDeg())-(g24VectTruth.Phi()*TMath::RadToDeg()));
          Bool_t Phi = kFALSE;

          if (diffPhi >= 170 && diffPhi <= 190) {
            Phi = kTRUE;		
	  }
		
	  if (Phi == kTRUE && Energy == kTRUE) {
	    ZDCHits = kTRUE;
	  }
        }
      }

    // Loop over our tracks, assign info
    for (SvtxTrackMap::Iter iter = trackmap->begin(); iter != trackmap->end();++iter)
      {
	SvtxTrack* track = iter->second;
	if ( track->get_pz() > 0 && track->get_charge() == 1){ // +ve z direction -> kaons, crappy way of selecting them for now w/o truth info
	  kaonVect.SetXYZ(track->get_px(), track->get_py(), track->get_pz());
	  kaonVectSmeared.SetXYZ(Position_Smear(track->get_px()), Position_Smear(track->get_py()), Position_Smear(track->get_pz()));
	  kaon4Vect.SetPxPyPzE(track->get_px(), track->get_py(), track->get_pz(), sqrt(pow(kaonVect.Mag(), 2)+pow(mKaonP,2)));
	  kaon4VectSmeared.SetPxPyPzE(kaonVectSmeared.x(), kaonVectSmeared.y(), kaonVectSmeared.z(), sqrt(pow(kaonVect.Mag(), 2)+pow(mKaonP,2)));
	}

	else if (track->get_pz() < 0  && track->get_charge() == -1 ){ // -ve z direction -> electrons, crappy way of selecting them for now w/o truth info
	  eVect.SetXYZ(track->get_px(), track->get_py(), track->get_pz());
	  eVectSmeared.SetXYZ(Position_Smear(track->get_px()), Position_Smear(track->get_py()), Position_Smear(track->get_pz()));
	  e4Vect.SetPxPyPzE(track->get_px(), track->get_py(), track->get_pz(), sqrt(pow(eVect.Mag(), 2)+pow(mElec,2)));
	  e4VectSmeared.SetPxPyPzE(eVectSmeared.x(), eVectSmeared.y(), eVectSmeared.z(), sqrt(pow(eVect.Mag(), 2)+pow(mElec,2)));
	}
      }

    // Now have relevant information from this event, fill some histograms and calculate some stuff

    // Calculate kinematic quantities
    virtphoton4Vect = eBeam4Vect - e4Vect;
    t4Vect = virtphoton4Vect - kaon4Vect;
    pmiss4Vect = (eBeam4Vect + pBeam4Vect) - (e4Vect+kaon4Vect);
    s4Vect = eBeam4Vect + pBeam4Vect; // added by Maggie Kerr, August 2, 2021 // are these the right 4 vectors?
    // s4Vect = virtphoton4Vect + pBeam4Vect; // added by Maggie Kerr, August 2, 2021 // are these the right 4 vectors?
    // s = s4Vect.Mag2();
    Q2 = -1*(virtphoton4Vect.Mag2());
    W = (virtphoton4Vect+pBeam4Vect).Mag();
    t = -(t4Vect.Mag2());
    xb =  Q2/(2*(pBeam4Vect.Dot(virtphoton4Vect)));
    xi = xb/(2-xb);
    // y = Q2/(xb*(s - mProt*mProt));
    // e = (2*(1-y))/(1 + (1-y)*(1-y));

    // Truth versions of kinematic quantities
    virtphoton4VectTruth = eBeam4Vect - e4VectTruth;
    t4VectTruth = virtphoton4VectTruth - kaon4VectTruth;
    pmiss4VectTruth = (eBeam4Vect + pBeam4Vect) - (e4VectTruth+kaon4VectTruth);
    Q2_truth = -1*(virtphoton4VectTruth.Mag2());
    W_truth = (virtphoton4VectTruth+pBeam4Vect).Mag();
    t_truth = -(t4VectTruth.Mag2());
    xb_truth =  Q2_truth/(2*(pBeam4Vect.Dot(virtphoton4VectTruth)));
    xi_truth = xb_truth/(2-xb_truth);

    // Fill Histograms

    ZDCHits = kTRUE;
    nHits = 3;
    
    if (ZDCHits == kTRUE) {

      h1_lTruth_E->Fill(l4VectTruth.E());
      h1_lTruth_P->Fill(l4VectTruth.P());
      h1_lTruth_Px->Fill(l4VectTruth.Px());
      h1_lTruth_Py->Fill(l4VectTruth.Py());
      h1_lTruth_Pz->Fill(l4VectTruth.Pz());
      h1_lTruth_Theta->Fill(l4VectTruth.Theta()*TMath::RadToDeg());
      h1_lTruth_Phi->Fill(l4VectTruth.Phi()*TMath::RadToDeg());
      h1_lTruth_PxPy->Fill(l4VectTruth.Px(),l4VectTruth.Py());
      h1_lTruth_ThetaP->Fill(l4VectTruth.Theta()*TMath::RadToDeg(),l4VectTruth.P());

      h1_eTruth_E->Fill(e4VectTruth.E());
      h1_eTruth_P->Fill(e4VectTruth.P());
      h1_eTruth_Px->Fill(e4VectTruth.Px());
      h1_eTruth_Py->Fill(e4VectTruth.Py());
      h1_eTruth_Pz->Fill(e4VectTruth.Pz());
      h1_eTruth_Theta->Fill(e4VectTruth.Theta()*TMath::RadToDeg());
      h1_eTruth_Phi->Fill(e4VectTruth.Phi()*TMath::RadToDeg());
      h1_eTruth_PxPy->Fill(e4VectTruth.Px(),e4VectTruth.Py());
      h1_eTruth_ThetaP->Fill(e4VectTruth.Theta()*TMath::RadToDeg(),e4VectTruth.P());

      h1_kTruth_E->Fill(kaon4VectTruth.E());
      h1_kTruth_P->Fill(kaon4VectTruth.P());
      h1_kTruth_Px->Fill(kaon4VectTruth.Px());
      h1_kTruth_Py->Fill(kaon4VectTruth.Py());
      h1_kTruth_Pz->Fill(kaon4VectTruth.Pz());
      h1_kTruth_Theta->Fill(kaon4VectTruth.Theta()*TMath::RadToDeg());
      h1_kTruth_Phi->Fill(kaon4VectTruth.Phi()*TMath::RadToDeg());
      h1_kTruth_PxPy->Fill(kaon4VectTruth.Px(),kaon4VectTruth.Py());
      h1_kTruth_ThetaP->Fill(kaon4VectTruth.Theta()*TMath::RadToDeg(),kaon4VectTruth.P());

      if (nHits == 3) {

        h1_nTruth_E->Fill(n4VectTruth.E());
        h1_nTruth_P->Fill(n4VectTruth.P());
        h1_nTruth_Px->Fill(n4VectTruth.Px());
        h1_nTruth_Py->Fill(n4VectTruth.Py());
        h1_nTruth_Pz->Fill(n4VectTruth.Pz());
        h1_nTruth_Theta->Fill(n4VectTruth.Theta()*TMath::RadToDeg());
        h1_nTruth_Phi->Fill(n4VectTruth.Phi()*TMath::RadToDeg());
        h1_nTruth_PxPy->Fill(n4VectTruth.Px(),n4VectTruth.Py());
        h1_nTruth_ThetaP->Fill(n4VectTruth.Theta()*TMath::RadToDeg(),n4VectTruth.P());

        h1_g1Truth_E->Fill(g14VectTruth.E());
        h1_g1Truth_P->Fill(g14VectTruth.P());
        h1_g1Truth_Px->Fill(g14VectTruth.Px());
        h1_g1Truth_Py->Fill(g14VectTruth.Py());
        h1_g1Truth_Pz->Fill(g14VectTruth.Pz());
        h1_g1Truth_Theta->Fill(g14VectTruth.Theta()*TMath::RadToDeg());
        h1_g1Truth_Phi->Fill(g14VectTruth.Phi()*TMath::RadToDeg());
        h1_g1Truth_PxPy->Fill(g14VectTruth.Px(),g14VectTruth.Py());
        h1_g1Truth_ThetaP->Fill(g14VectTruth.Theta()*TMath::RadToDeg(),g14VectTruth.P());
    
        h1_g2Truth_E->Fill(g24VectTruth.E());
        h1_g2Truth_P->Fill(g24VectTruth.P());
        h1_g2Truth_Px->Fill(g24VectTruth.Px());
        h1_g2Truth_Py->Fill(g24VectTruth.Py());
        h1_g2Truth_Pz->Fill(g24VectTruth.Pz());
        h1_g2Truth_Theta->Fill(g24VectTruth.Theta()*TMath::RadToDeg());
        h1_g2Truth_Phi->Fill(g24VectTruth.Phi()*TMath::RadToDeg());
        h1_g2Truth_PxPy->Fill(g24VectTruth.Px(),g24VectTruth.Py());
        h1_g2Truth_ThetaP->Fill(g24VectTruth.Theta()*TMath::RadToDeg(),g24VectTruth.P());
      }

      h1_Q2_Dist->Fill(Q2);
      h1_W_Dist->Fill(W);
      h1_t_Dist->Fill(t);
      h1_xb_Dist->Fill(xb);
      h1_xi_Dist->Fill(xi);

      h1_Q2Truth_Dist->Fill(((Q2-Q2_truth)/Q2_truth)*100);
      h1_WTruth_Dist->Fill(((W-W_truth)/W_truth)*100);
      h1_tTruth_Dist->Fill(((t-t_truth)/t_truth)*100);
      h1_xbTruth_Dist->Fill(((xb-xb_truth)/xb_truth)*100);
      h1_xiTruth_Dist->Fill(((xi-xi_truth)/xi_truth)*100);

      if (Q2 >= 0 && Q2 < 2.5) {   
       	h2_tTruthvQ2_0_25_Dist->Fill(((t-t_truth)/t_truth)*100, Q2);
      }

      else if (Q2 >= 2.5 && Q2 < 3.75) {   
       	h2_tTruthvQ2_25_375_Dist->Fill(((t-t_truth)/t_truth)*100, Q2);
      }

      else if (Q2 >= 3.75 && Q2 < 5) {   
       	h2_tTruthvQ2_375_5_Dist->Fill(((t-t_truth)/t_truth)*100, Q2);
      }

      else if (Q2 >= 5 && Q2 < 6.25) {   
       	h2_tTruthvQ2_5_625_Dist->Fill(((t-t_truth)/t_truth)*100, Q2);
      }

      else if (Q2 >= 6.25 && Q2 < 7.5) {   
       	h2_tTruthvQ2_625_75_Dist->Fill(((t-t_truth)/t_truth)*100, Q2);
      }

      else if (Q2 >= 7.5 && Q2 < 8.75) {   
       	h2_tTruthvQ2_75_875_Dist->Fill(((t-t_truth)/t_truth)*100, Q2);
      }

      else if (Q2 >= 8.75 && Q2 < 10) {   
       	h2_tTruthvQ2_875_10_Dist->Fill(((t-t_truth)/t_truth)*100, Q2);
      }

      else if (Q2 >= 10 && Q2 < 11.25) {   
       	h2_tTruthvQ2_10_1125_Dist->Fill(((t-t_truth)/t_truth)*100, Q2);
      }

      else if (Q2 >= 11.25) {   
       	h2_tTruthvQ2_g1125_Dist->Fill(((t-t_truth)/t_truth)*100, Q2);
      }

    h1_k_E->Fill(kaon4Vect.E());
    h1_k_P->Fill(kaon4Vect.P());
    h1_k_Px->Fill(kaon4Vect.Px());
    h1_k_Py->Fill(kaon4Vect.Py());
    h1_k_Pz->Fill(kaon4Vect.Pz());
    h1_k_Theta->Fill(kaon4Vect.Theta()*TMath::RadToDeg());
    h1_k_Phi->Fill(kaon4Vect.Phi()*TMath::RadToDeg());
    h2_k_PxPy->Fill(kaon4Vect.Px(),kaon4Vect.Py());
    h2_k_ThetaP->Fill(kaon4Vect.Theta()*TMath::RadToDeg(),kaon4Vect.P());

    h1_e_E->Fill(e4Vect.E());
    h1_e_P->Fill(e4Vect.P());
    h1_e_Px->Fill(e4Vect.Px());
    h1_e_Py->Fill(e4Vect.Py());
    h1_e_Pz->Fill(e4Vect.Pz());
    h1_e_Theta->Fill(e4Vect.Theta()*TMath::RadToDeg());
    h1_e_Phi->Fill(e4Vect.Phi()*TMath::RadToDeg());
    h2_e_PxPy->Fill(e4Vect.Px(),e4Vect.Py());
    h2_e_ThetaP->Fill(e4Vect.Theta()*TMath::RadToDeg(),e4Vect.P());

    h1_kTruthC_E->Fill(((kaon4Vect.E()-kaon4VectTruth.E())/kaon4VectTruth.E())*100);
    h1_kTruthC_P->Fill(((kaon4Vect.P()-kaon4VectTruth.P())/kaon4VectTruth.P())*100);
    h1_kTruthC_Px->Fill(((kaon4Vect.Px()-kaon4VectTruth.Px())/kaon4VectTruth.Px())*100);
    h1_kTruthC_Py->Fill(((kaon4Vect.Py()-kaon4VectTruth.Py())/kaon4VectTruth.Py())*100);
    h1_kTruthC_Pz->Fill(((kaon4Vect.Pz()-kaon4VectTruth.Pz())/kaon4VectTruth.Pz())*100);
    h1_kTruthC_Theta->Fill(((kaon4Vect.Theta()*TMath::RadToDeg()-kaon4VectTruth.Theta()*TMath::RadToDeg())/kaon4VectTruth.Theta()*TMath::RadToDeg())*100);
    h1_kTruthC_Phi->Fill(((kaon4Vect.Phi()*TMath::RadToDeg()-kaon4VectTruth.Phi()*TMath::RadToDeg())/kaon4VectTruth.Phi()*TMath::RadToDeg())*100);
    h2_kTruthC_PxPy->Fill(((kaon4Vect.Px()-kaon4VectTruth.Px())/kaon4VectTruth.Px())*100,((kaon4Vect.Py()-kaon4VectTruth.Py())/kaon4VectTruth.Py())*100);
    h2_kTruthC_ThetaP->Fill(((kaon4Vect.Theta()*TMath::RadToDeg()-kaon4VectTruth.Theta()*TMath::RadToDeg())/kaon4VectTruth.Theta()*TMath::RadToDeg())*100,((kaon4Vect.P()-kaon4VectTruth.P())/kaon4VectTruth.P())*100);

    h1_eTruthC_E->Fill(((e4Vect.E()-e4VectTruth.E())/e4VectTruth.E())*100);
    h1_eTruthC_P->Fill(((e4Vect.P()-e4VectTruth.P())/e4VectTruth.P())*100);
    h1_eTruthC_Px->Fill(((e4Vect.Px()-e4VectTruth.Px())/e4VectTruth.Px())*100);
    h1_eTruthC_Py->Fill(((e4Vect.Py()-e4VectTruth.Py())/e4VectTruth.Py())*100);
    h1_eTruthC_Pz->Fill(((e4Vect.Pz()-e4VectTruth.Pz())/e4VectTruth.Pz())*100);
    h1_eTruthC_Theta->Fill(((e4Vect.Theta()*TMath::RadToDeg()-e4VectTruth.Theta()*TMath::RadToDeg())/e4VectTruth.Theta()*TMath::RadToDeg())*100);
    h1_eTruthC_Phi->Fill(((e4Vect.Phi()*TMath::RadToDeg()-e4VectTruth.Phi()*TMath::RadToDeg())/e4VectTruth.Phi()*TMath::RadToDeg())*100);
    h2_eTruthC_PxPy->Fill(((e4Vect.Px()-e4VectTruth.Px())/e4VectTruth.Px())*100,((e4Vect.Py()-e4VectTruth.Py())/e4VectTruth.Py())*100);
    h2_eTruthC_ThetaP->Fill(((e4Vect.Theta()*TMath::RadToDeg()-e4VectTruth.Theta()*TMath::RadToDeg())/e4VectTruth.Theta()*TMath::RadToDeg())*100,((e4Vect.P()-e4VectTruth.P())/e4VectTruth.P())*100);

    h2_ZDC_XY->Fill(nZDCPos.x(), nZDCPos.y());
    h2_ZDC_XY_Smeared->Fill(nZDCPosSmeared.x(), nZDCPosSmeared.y());


    // Fill weighted histograms
    h1_K_px_Weighted->Fill(kaon4Vect.Px(), wgt);
    h1_K_py_Weighted->Fill(kaon4Vect.Py(), wgt);
    h1_K_pz_Weighted->Fill(kaon4Vect.Pz(), wgt);
    h1_K_p_Weighted->Fill(kaon4Vect.P(), wgt);
    h1_K_E_Weighted->Fill(kaon4Vect.E(), wgt);
    h1_K_Theta_Weighted->Fill(kaon4Vect.Theta()*TMath::RadToDeg(), wgt);
    h1_K_Phi_Weighted->Fill(kaon4Vect.Phi()*TMath::RadToDeg(), wgt);
    h1_e_px_Weighted->Fill(e4Vect.Px(), wgt);
    h1_e_py_Weighted->Fill(e4Vect.Py(), wgt);
    h1_e_pz_Weighted->Fill(e4Vect.Pz(), wgt);
    h1_e_p_Weighted->Fill(e4Vect.P(), wgt);
    h1_e_E_Weighted->Fill(e4Vect.E(), wgt);
    h1_e_Theta_Weighted->Fill(e4Vect.Theta()*TMath::RadToDeg(), wgt);
    h1_e_Phi_Weighted->Fill(e4Vect.Phi()*TMath::RadToDeg(), wgt);
    h1_L_px_Weighted->Fill(l4Vect.Px(), wgt);
    h1_L_py_Weighted->Fill(l4Vect.Py(), wgt);
    h1_L_pz_Weighted->Fill(l4Vect.Pz(), wgt);
    h1_L_p_Weighted->Fill(l4Vect.P(), wgt);
    h1_L_E_Weighted->Fill(l4Vect.E(), wgt);
    h1_L_Theta_Weighted->Fill(l4Vect.Theta()*TMath::RadToDeg(), wgt);
    h1_L_Phi_Weighted->Fill(l4Vect.Phi()*TMath::RadToDeg(), wgt);

    if (nHits == 3) {

      h1_n_px_Weighted->Fill(n4Vect.Px(), wgt);
      h1_n_py_Weighted->Fill(n4Vect.Py(), wgt);
      h1_n_pz_Weighted->Fill(n4Vect.Pz(), wgt);
      h1_n_p_Weighted->Fill(n4Vect.P(), wgt);
      h1_n_E_Weighted->Fill(n4Vect.E(), wgt);
      h1_n_Theta_Weighted->Fill(n4Vect.Theta()*TMath::RadToDeg(), wgt);
      h1_n_Phi_Weighted->Fill(n4Vect.Phi()*TMath::RadToDeg(), wgt);

      h1_g1_px_Weighted->Fill(g14Vect.Px(), wgt);
      h1_g1_py_Weighted->Fill(g14Vect.Py(), wgt);
      h1_g1_pz_Weighted->Fill(g14Vect.Pz(), wgt);
      h1_g1_p_Weighted->Fill(g14Vect.P(), wgt);
      h1_g1_E_Weighted->Fill(g14Vect.E(), wgt);
      h1_g1_Theta_Weighted->Fill(g14Vect.Theta()*TMath::RadToDeg(), wgt);   
      h1_g1_Phi_Weighted->Fill(g14Vect.Phi()*TMath::RadToDeg(), wgt);

      h1_g2_px_Weighted->Fill(g24Vect.Px(), wgt);
      h1_g2_py_Weighted->Fill(g24Vect.Py(), wgt);
      h1_g2_pz_Weighted->Fill(g24Vect.Pz(), wgt);
      h1_g2_p_Weighted->Fill(g24Vect.P(), wgt);
      h1_g2_E_Weighted->Fill(g24Vect.E(), wgt);
      h1_g2_Theta_Weighted->Fill(g24Vect.Theta()*TMath::RadToDeg(), wgt);
      h1_g2_Phi_Weighted->Fill(g24Vect.Phi()*TMath::RadToDeg(), wgt);

    }

    h1_Q2_Dist_Weighted->Fill(Q2, wgt);
    h1_W_Dist_Weighted->Fill(W, wgt);
    h1_t_Dist_Weighted->Fill(t, wgt);
    h1_t_alt_Dist_Weighted->Fill(t_alt, wgt);
    h1_t_comp_Weighted->Fill(((t_alt-t)/t)*100, wgt);
    h1_xb_Dist_Weighted->Fill(xb, wgt);
    h1_xi_Dist_Weighted->Fill(xi, wgt);

    h1_Q2Truth_Dist_Weighted->Fill(Q2_truth, wgt);
    h1_WTruth_Dist_Weighted->Fill(W_truth, wgt);
    h1_tTruth_Dist_Weighted->Fill(t_truth, wgt);
    h1_xbTruth_Dist_Weighted->Fill(xb_truth, wgt);
    h1_xiTruth_Dist_Weighted->Fill(xi_truth, wgt);

    h2_t_ep_Weighted->Fill(t, e4Vect.P(), wgt);
    h2_t_Q2_Weighted->Fill(t,Q2, wgt);
    h2_delta_t_t_Weighted->Fill(((t - t_truth)/t_truth)*100, t, wgt);
    
    for(Int_t B = 0; B < 7; B++){
      Double_t Q2_low = 5+(B*5);
      Double_t Q2_high = 10+(B*5);
      if ( Q2_truth > Q2_low && Q2_truth < Q2_high){
	h1_t_Q2_Weighted[B]->Fill(t, wgt);
	h1_t_alt_Q2_Weighted[B]->Fill(t_alt, wgt);
	h2_delta_t_t_Q2_Weighted[B]->Fill(((t - t_truth)/t_truth)*100, t, wgt);
      }

      h1_KTruth_p_Weighted->Fill((kaon4Vect.P()-kaon4VectTruth.P())/(kaon4VectTruth.P())*100, wgt);
      h1_KTruth_px_Weighted->Fill((kaon4Vect.Px()-kaon4VectTruth.Px())/(kaon4VectTruth.Px())*100, wgt);
      h1_KTruth_py_Weighted->Fill((kaon4Vect.Py()-kaon4VectTruth.Py())/(kaon4VectTruth.Py())*100, wgt);
      h1_KTruth_pz_Weighted->Fill((kaon4Vect.Pz()-kaon4VectTruth.Pz())/(kaon4VectTruth.Pz())*100, wgt);
      h1_KTruth_E_Weighted->Fill((kaon4Vect.E()-kaon4VectTruth.E())/(kaon4VectTruth.E())*100, wgt);
      h1_eTruth_p_Weighted->Fill((e4Vect.P()-e4VectTruth.P())/(e4VectTruth.P())*100, wgt);
      h1_eTruth_px_Weighted->Fill((e4Vect.Px()-e4VectTruth.Px())/(e4VectTruth.Px())*100, wgt);
      h1_eTruth_py_Weighted->Fill((e4Vect.Py()-e4VectTruth.Py())/(e4VectTruth.Py())*100, wgt);
      h1_eTruth_pz_Weighted->Fill((e4Vect.Pz()-e4VectTruth.Pz())/(e4VectTruth.Pz())*100, wgt);
      h1_eTruth_E_Weighted->Fill((e4Vect.E()-e4VectTruth.E())/(e4VectTruth.E())*100, wgt);
      h1_LTruth_p_Weighted->Fill((l4Vect.P()-l4VectTruth.P())/(l4VectTruth.P())*100, wgt);
      h1_LTruth_px_Weighted->Fill((l4Vect.Px()-l4VectTruth.Px())/(l4VectTruth.Px())*100, wgt);
      h1_LTruth_py_Weighted->Fill((l4Vect.Py()-l4VectTruth.Py())/(l4VectTruth.Py())*100, wgt);
      h1_LTruth_pz_Weighted->Fill((l4Vect.Pz()-l4VectTruth.Pz())/(l4VectTruth.Pz())*100, wgt);
      h1_LTruth_E_Weighted->Fill((l4Vect.E()-l4VectTruth.E())/(l4VectTruth.E())*100, wgt);

      if (nHits == 3) {

	h1_nTruth_p_Weighted->Fill((n4Vect.P()-n4VectTruth.P())/(n4VectTruth.P())*100, wgt);
	h1_nTruth_px_Weighted->Fill((n4Vect.Px()-n4VectTruth.Px())/(n4VectTruth.Px())*100, wgt);
	h1_nTruth_py_Weighted->Fill((n4Vect.Py()-n4VectTruth.Py())/(n4VectTruth.Py())*100, wgt);
	h1_nTruth_pz_Weighted->Fill((n4Vect.Pz()-n4VectTruth.Pz())/(n4VectTruth.Pz())*100, wgt);
	h1_nTruth_E_Weighted->Fill((n4Vect.E()-n4VectTruth.E())/(n4VectTruth.E())*100, wgt);

	h1_g1Truth_p_Weighted->Fill((g14Vect.P()-g14VectTruth.P())/(g14VectTruth.P())*100, wgt);
	h1_g1Truth_px_Weighted->Fill((g14Vect.Px()-g14VectTruth.Px())/(g14VectTruth.Px())*100, wgt);
	h1_g1Truth_py_Weighted->Fill((g14Vect.Py()-g14VectTruth.Py())/(g14VectTruth.Py())*100, wgt);
	h1_g1Truth_pz_Weighted->Fill((g14Vect.Pz()-g14VectTruth.Pz())/(g14VectTruth.Pz())*100, wgt);
	h1_g1Truth_E_Weighted->Fill((g14Vect.E()-g14VectTruth.E())/(g14VectTruth.E())*100, wgt);

	h1_g2Truth_p_Weighted->Fill((g24Vect.P()-g24VectTruth.P())/(g24VectTruth.P())*100, wgt);
	h1_g2Truth_px_Weighted->Fill((g24Vect.Px()-g24VectTruth.Px())/(g24VectTruth.Px())*100, wgt);
	h1_g2Truth_py_Weighted->Fill((g24Vect.Py()-g24VectTruth.Py())/(g24VectTruth.Py())*100, wgt);
	h1_g2Truth_pz_Weighted->Fill((g24Vect.Pz()-g24VectTruth.Pz())/(g24VectTruth.Pz())*100, wgt);
	h1_g2Truth_E_Weighted->Fill((g24Vect.E()-g24VectTruth.E())/(g24VectTruth.E())*100, wgt);

      }

      h1_KTruth_p_Smeared_Weighted->Fill((kaon4VectSmeared.P()-kaon4VectTruth.P())/(kaon4VectTruth.P())*100, wgt);
      h1_KTruth_px_Smeared_Weighted->Fill((kaon4VectSmeared.Px()-kaon4VectTruth.Px())/(kaon4VectTruth.Px())*100, wgt);
      h1_KTruth_py_Smeared_Weighted->Fill((kaon4VectSmeared.Py()-kaon4VectTruth.Py())/(kaon4VectTruth.Py())*100, wgt);
      h1_KTruth_pz_Smeared_Weighted->Fill((kaon4VectSmeared.Pz()-kaon4VectTruth.Pz())/(kaon4VectTruth.Pz())*100, wgt);
      h1_KTruth_E_Smeared_Weighted->Fill((kaon4VectSmeared.E()-kaon4VectTruth.E())/(kaon4VectTruth.E())*100, wgt);
      h1_eTruth_p_Smeared_Weighted->Fill((e4VectSmeared.P()-e4VectTruth.P())/(e4VectTruth.P())*100, wgt);
      h1_eTruth_px_Smeared_Weighted->Fill((e4VectSmeared.Px()-e4VectTruth.Px())/(e4VectTruth.Px())*100, wgt);
      h1_eTruth_py_Smeared_Weighted->Fill((e4VectSmeared.Py()-e4VectTruth.Py())/(e4VectTruth.Py())*100, wgt);
      h1_eTruth_pz_Smeared_Weighted->Fill((e4VectSmeared.Pz()-e4VectTruth.Pz())/(e4VectTruth.Pz())*100, wgt);
      h1_eTruth_E_Smeared_Weighted->Fill((e4VectSmeared.E()-e4VectTruth.E())/(e4VectTruth.E())*100, wgt);
      h1_LTruth_p_Smeared_Weighted->Fill((l4VectSmeared.P()-l4VectTruth.P())/(l4VectTruth.P())*100, wgt);
      h1_LTruth_px_Smeared_Weighted->Fill((l4VectSmeared.Px()-l4VectTruth.Px())/(l4VectTruth.Px())*100, wgt);
      h1_LTruth_py_Smeared_Weighted->Fill((l4VectSmeared.Py()-l4VectTruth.Py())/(l4VectTruth.Py())*100, wgt);
      h1_LTruth_pz_Smeared_Weighted->Fill((l4VectSmeared.Pz()-l4VectTruth.Pz())/(l4VectTruth.Pz())*100, wgt);
      h1_LTruth_E_Smeared_Weighted->Fill((l4VectSmeared.E()-l4VectTruth.E())/(l4VectTruth.E())*100, wgt);

      if (nHits == 3) {

	h1_nTruth_p_Smeared_Weighted->Fill((n4VectSmeared.P()-n4VectTruth.P())/(n4VectTruth.P())*100, wgt);
	h1_nTruth_px_Smeared_Weighted->Fill((n4VectSmeared.Px()-n4VectTruth.Px())/(n4VectTruth.Px())*100, wgt);
	h1_nTruth_py_Smeared_Weighted->Fill((n4VectSmeared.Py()-n4VectTruth.Py())/(n4VectTruth.Py())*100, wgt);
	h1_nTruth_pz_Smeared_Weighted->Fill((n4VectSmeared.Pz()-n4VectTruth.Pz())/(n4VectTruth.Pz())*100, wgt);
	h1_nTruth_E_Smeared_Weighted->Fill((n4VectSmeared.E()-n4VectTruth.E())/(n4VectTruth.E())*100, wgt);

	h1_g1Truth_p_Smeared_Weighted->Fill((g14VectSmeared.P()-g14VectTruth.P())/(g14VectTruth.P())*100, wgt);
	h1_g1Truth_px_Smeared_Weighted->Fill((g14VectSmeared.Px()-g14VectTruth.Px())/(g14VectTruth.Px())*100, wgt);
	h1_g1Truth_py_Smeared_Weighted->Fill((g14VectSmeared.Py()-g14VectTruth.Py())/(g14VectTruth.Py())*100, wgt);
	h1_g1Truth_pz_Smeared_Weighted->Fill((g14VectSmeared.Pz()-g14VectTruth.Pz())/(g14VectTruth.Pz())*100, wgt);
	h1_g1Truth_E_Smeared_Weighted->Fill((g14VectSmeared.E()-g14VectTruth.E())/(g14VectTruth.E())*100, wgt);

	h1_g2Truth_p_Smeared_Weighted->Fill((g24VectSmeared.P()-g24VectTruth.P())/(g24VectTruth.P())*100, wgt);
	h1_g2Truth_px_Smeared_Weighted->Fill((g24VectSmeared.Px()-g24VectTruth.Px())/(g24VectTruth.Px())*100, wgt);
	h1_g2Truth_py_Smeared_Weighted->Fill((g24VectSmeared.Py()-g24VectTruth.Py())/(g24VectTruth.Py())*100, wgt);
	h1_g2Truth_pz_Smeared_Weighted->Fill((g24VectSmeared.Pz()-g24VectTruth.Pz())/(g24VectTruth.Pz())*100, wgt);
	h1_g2Truth_E_Smeared_Weighted->Fill((g24VectSmeared.E()-g24VectTruth.E())/(g24VectTruth.E())*100, wgt);

      }
    
      h2_ZDC_XY_Weighted->Fill(nZDCPos.x(), nZDCPos.y(), wgt);
      h2_ZDC_XY_Smeared_Weighted->Fill(nZDCPosSmeared.x(), nZDCPosSmeared.y(), wgt);

      h2_KTrack_ThetaPhi_Weighted->Fill((kaon4Vect.Theta()*TMath::RadToDeg()), (kaon4Vect.Phi()*TMath::RadToDeg()), wgt);
      h2_KTrack_pTheta_Weighted->Fill((kaon4Vect.Theta()*TMath::RadToDeg()), kaon4Vect.P(), wgt);
      h2_eTrack_ThetaPhi_Weighted->Fill((e4Vect.Theta()*TMath::RadToDeg()), (e4Vect.Phi()*TMath::RadToDeg()), wgt);
      h2_eTrack_pTheta_Weighted->Fill((e4Vect.Theta()*TMath::RadToDeg()), e4Vect.P(), wgt);
      h2_LTrack_ThetaPhi_Weighted->Fill((l4Vect.Theta()*TMath::RadToDeg()), (l4Vect.Phi()*TMath::RadToDeg()), wgt);
      h2_LTrack_pTheta_Weighted->Fill((l4Vect.Theta()*TMath::RadToDeg()), l4Vect.P(), wgt);
      h2_KTrack_ThetaPhi_Smeared_Weighted->Fill((kaon4VectSmeared.Theta()*TMath::RadToDeg()), (kaon4VectSmeared.Phi()*TMath::RadToDeg()), wgt);
      h2_KTrack_pTheta_Smeared_Weighted->Fill((kaon4VectSmeared.Theta()*TMath::RadToDeg()), kaon4VectSmeared.P(), wgt);
      h2_eTrack_ThetaPhi_Smeared_Weighted->Fill((e4VectSmeared.Theta()*TMath::RadToDeg()), (e4VectSmeared.Phi()*TMath::RadToDeg()), wgt);
      h2_eTrack_pTheta_Smeared_Weighted->Fill((e4VectSmeared.Theta()*TMath::RadToDeg()), e4VectSmeared.P(), wgt);
      h2_LTrack_ThetaPhi_Smeared_Weighted->Fill((l4VectSmeared.Theta()*TMath::RadToDeg()), (l4VectSmeared.Phi()*TMath::RadToDeg()), wgt);
      h2_LTrack_pTheta_Smeared_Weighted->Fill((l4VectSmeared.Theta()*TMath::RadToDeg()), l4VectSmeared.P(), wgt);

      if (nHits == 3) {

	h2_nTrack_ThetaPhi_Weighted->Fill((n4Vect.Theta()*TMath::RadToDeg()), (n4Vect.Phi()*TMath::RadToDeg()), wgt);
	h2_nTrack_pTheta_Weighted->Fill((n4Vect.Theta()*TMath::RadToDeg()), n4Vect.P(), wgt);

	h2_g1Track_ThetaPhi_Weighted->Fill((g14Vect.Theta()*TMath::RadToDeg()), (g14Vect.Phi()*TMath::RadToDeg()), wgt);
	h2_g1Track_pTheta_Weighted->Fill((g14Vect.Theta()*TMath::RadToDeg()), g14Vect.P(), wgt);

	h2_g2Track_ThetaPhi_Weighted->Fill((g24Vect.Theta()*TMath::RadToDeg()), (g24Vect.Phi()*TMath::RadToDeg()), wgt);
	h2_g2Track_pTheta_Weighted->Fill((g24Vect.Theta()*TMath::RadToDeg()), g24Vect.P(), wgt);

	h2_nTrack_ThetaPhi_Smeared_Weighted->Fill((n4VectSmeared.Theta()*TMath::RadToDeg()), (n4VectSmeared.Phi()*TMath::RadToDeg()), wgt);
	h2_nTrack_pTheta_Smeared_Weighted->Fill((n4VectSmeared.Theta()*TMath::RadToDeg()), n4VectSmeared.P(), wgt);

	h2_g1Track_ThetaPhi_Smeared_Weighted->Fill((g14VectSmeared.Theta()*TMath::RadToDeg()), (g14VectSmeared.Phi()*TMath::RadToDeg()), wgt);
	h2_g1Track_pTheta_Smeared_Weighted->Fill((g14VectSmeared.Theta()*TMath::RadToDeg()), g14VectSmeared.P(), wgt);

	h2_g2Track_ThetaPhi_Smeared_Weighted->Fill((g24VectSmeared.Theta()*TMath::RadToDeg()), (g24VectSmeared.Phi()*TMath::RadToDeg()), wgt);
	h2_g2Track_pTheta_Smeared_Weighted->Fill((g24VectSmeared.Theta()*TMath::RadToDeg()), g24VectSmeared.P(), wgt);

      }
    
      h2_KTruth_pxpy_Weighted->Fill((kaon4Vect.Px()-kaon4VectTruth.Px())/(kaon4VectTruth.Px())*100, (kaon4Vect.Py()-kaon4VectTruth.Py())/(kaon4VectTruth.Py())*100, wgt);
      h2_eTruth_pxpy_Weighted->Fill((e4Vect.Px()-e4VectTruth.Px())/(e4VectTruth.Px())*100, (e4Vect.Py()-e4VectTruth.Py())/(e4VectTruth.Py())*100, wgt);
      h2_LTruth_pxpy_Weighted->Fill((l4Vect.Px()-l4VectTruth.Px())/(l4VectTruth.Px())*100, (l4Vect.Py()-l4VectTruth.Py())/(l4VectTruth.Py())*100, wgt);

      h2_KTruth_pxpz_Weighted->Fill((kaon4Vect.Px()-kaon4VectTruth.Px())/(kaon4VectTruth.Px())*100, (kaon4Vect.Pz()-kaon4VectTruth.Pz())/(kaon4VectTruth.Pz())*100, wgt);
      h2_eTruth_pxpz_Weighted->Fill((e4Vect.Px()-e4VectTruth.Px())/(e4VectTruth.Px())*100, (e4Vect.Pz()-e4VectTruth.Pz())/(e4VectTruth.Pz())*100, wgt);
      h2_LTruth_pxpz_Weighted->Fill((l4Vect.Px()-l4VectTruth.Px())/(l4VectTruth.Px())*100, (l4Vect.Pz()-l4VectTruth.Pz())/(l4VectTruth.Pz())*100, wgt);

      h2_KTruth_pypz_Weighted->Fill((kaon4Vect.Py()-kaon4VectTruth.Py())/(kaon4VectTruth.Py())*100, (kaon4Vect.Pz()-kaon4VectTruth.Pz())/(kaon4VectTruth.Pz())*100, wgt);
      h2_eTruth_pypz_Weighted->Fill((e4Vect.Py()-e4VectTruth.Py())/(e4VectTruth.Py())*100, (e4Vect.Pz()-e4VectTruth.Pz())/(e4VectTruth.Pz())*100, wgt);
      h2_LTruth_pypz_Weighted->Fill((l4Vect.Py()-l4VectTruth.Py())/(l4VectTruth.Py())*100, (l4Vect.Pz()-l4VectTruth.Pz())/(l4VectTruth.Pz())*100, wgt);

      h2_KTruth_pxpy_Smeared_Weighted->Fill((kaon4VectSmeared.Px()-kaon4VectTruth.Px())/(kaon4VectTruth.Px())*100, (kaon4VectSmeared.Py()-kaon4VectTruth.Py())/(kaon4VectTruth.Py())*100, wgt);
      h2_eTruth_pxpy_Smeared_Weighted->Fill((e4VectSmeared.Px()-e4VectTruth.Px())/(e4VectTruth.Px())*100, (e4VectSmeared.Py()-e4VectTruth.Py())/(e4VectTruth.Py())*100, wgt);
      h2_LTruth_pxpy_Smeared_Weighted->Fill((l4VectSmeared.Px()-l4VectTruth.Px())/(l4VectTruth.Px())*100, (l4VectSmeared.Py()-l4VectTruth.Py())/(l4VectTruth.Py())*100, wgt);

      if (nHits == 3) {

	h2_nTruth_pxpy_Weighted->Fill((n4Vect.Px()-n4VectTruth.Px())/(n4VectTruth.Px())*100, (n4Vect.Py()-n4VectTruth.Py())/(n4VectTruth.Py())*100, wgt);
	h2_g1Truth_pxpy_Weighted->Fill((g14Vect.Px()-g14VectTruth.Px())/(g14VectTruth.Px())*100, (g14Vect.Py()-g14VectTruth.Py())/(g14VectTruth.Py())*100, wgt);
	h2_g2Truth_pxpy_Weighted->Fill((g24Vect.Px()-g24VectTruth.Px())/(g24VectTruth.Px())*100, (g24Vect.Py()-g24VectTruth.Py())/(g24VectTruth.Py())*100, wgt);

	h2_nTruth_pxpz_Weighted->Fill((n4Vect.Px()-n4VectTruth.Px())/(n4VectTruth.Px())*100, (n4Vect.Pz()-n4VectTruth.Pz())/(n4VectTruth.Pz())*100, wgt);
	h2_g1Truth_pxpz_Weighted->Fill((g14Vect.Px()-g14VectTruth.Px())/(g14VectTruth.Px())*100, (g14Vect.Pz()-g14VectTruth.Pz())/(g14VectTruth.Pz())*100, wgt);
	h2_g2Truth_pxpz_Weighted->Fill((g24Vect.Px()-g24VectTruth.Px())/(g24VectTruth.Px())*100, (g24Vect.Pz()-g24VectTruth.Pz())/(g24VectTruth.Pz())*100, wgt);

	h2_nTruth_pypz_Weighted->Fill((n4Vect.Py()-n4VectTruth.Py())/(n4VectTruth.Py())*100, (n4Vect.Pz()-n4VectTruth.Pz())/(n4VectTruth.Pz())*100, wgt);
	h2_g1Truth_pypz_Weighted->Fill((g14Vect.Py()-g14VectTruth.Py())/(g14VectTruth.Py())*100, (g14Vect.Pz()-g14VectTruth.Pz())/(g14VectTruth.Pz())*100, wgt);
	h2_g2Truth_pypz_Weighted->Fill((g24Vect.Py()-g24VectTruth.Py())/(g24VectTruth.Py())*100, (g24Vect.Pz()-g24VectTruth.Pz())/(g24VectTruth.Pz())*100, wgt);

	h2_nTruth_pxpy_Smeared_Weighted->Fill((n4VectSmeared.Px()-n4VectTruth.Px())/(n4VectTruth.Px())*100, (n4VectSmeared.Py()-n4VectTruth.Py())/(n4VectTruth.Py())*100, wgt);
	h2_g1Truth_pxpy_Smeared_Weighted->Fill((g14VectSmeared.Px()-g14VectTruth.Px())/(g14VectTruth.Px())*100, (g14VectSmeared.Py()-g14VectTruth.Py())/(g14VectTruth.Py())*100, wgt);
	h2_g2Truth_pxpy_Smeared_Weighted->Fill((g24VectSmeared.Px()-g24VectTruth.Px())/(g24VectTruth.Px())*100, (g24VectSmeared.Py()-g24VectTruth.Py())/(g24VectTruth.Py())*100, wgt);

      }

    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
  
//____________________________________________________________________________..
int ECCE_kLambda::ResetEvent(PHCompositeNode *topNode)
{
//  std::cout << "ECCE_kLambda::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
//
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_kLambda::EndRun(const int runnumber)
{
  std::cout << "ECCE_kLambda::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_kLambda::End(PHCompositeNode *topNode)
{
  std::cout << "ECCE_kLambda::End(PHCompositeNode *topNode) This is the End..." << std::endl;

  outfile->cd();
  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(outfilename, "UPDATE");

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_kLambda::Reset(PHCompositeNode *topNode)
{
 std::cout << "ECCE_kLambda::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void ECCE_kLambda::Print(const std::string &what) const
{
  std::cout << "ECCE_kLambda::Print(const std::string &what) const Printing info for " << what << std::endl;
}

//***************************************************

float ECCE_kLambda::EMCAL_Smear(float E) {

  float resolution, E_reco;

  resolution = sqrt(.45*.45/E + 0.075*0.075);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;
}

//*****************************************************

float ECCE_kLambda::HCAL_Smear(float E) {

  float resolution, E_reco;

  resolution = sqrt(.50*.50/E + 0.1*0.1);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;
}

//*****************************************************

float ECCE_kLambda::PbWO4_Smear(float E) {

  float resolution, E_reco;

  resolution = sqrt(.25*.25/E + 0.04*0.04);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;

}

//*****************************************************

float ECCE_kLambda::Position_Smear(float P) {

  float resolution, P_reco;

  resolution = 0.1;         /// Position resolution 0.1 cm
  P_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * P;

  return P_reco;

}

//*****************************************************

bool ECCE_kLambda::Check_eKaon(PHCompositeNode* topNode)
{
  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!trackmap)
    {
      trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
      if (!trackmap)
    	{
    	  cout << "ECCE_kLambda::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
    	  exit(-1);
    	}
    }
  int nTracks = 0;
  Bool_t ElecTrack = kFALSE;
  Bool_t KaonTrack = kFALSE;
  // Iterate over tracks
  for (SvtxTrackMap::Iter iter = trackmap->begin(); iter != trackmap->end(); ++iter){
    SvtxTrack* track = iter->second;
    nTracks++;
    // std::cout << track->get_charge() << endl;
    if ( track->get_pz() > 0 && track->get_charge() == 1){ // +ve z direction -> kaons, crappy way of selecting them for now w/o truth info
      KaonTrack = kTRUE;
    }
    else if (track->get_pz() < 0  && track->get_charge() == -1 ){ // -ve z direction -> electrons, crappy way of selecting them for now w/o truth info
      ElecTrack = kTRUE;
     }
    }
  
  // std::cout << nTracks << endl;
  if( KaonTrack == kTRUE && ElecTrack == kTRUE && nTracks == 2){ // Both a kaon and an electron track, only 2 tracks
    return true;
  }
  else{
    return false;
  }
}


//*****************************************************

bool ECCE_kLambda::Check_hits(PHCompositeNode* topNode)
{
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
  if (!truthinfo) { 
    cout << PHWHERE << "PHG4TruthInfoContainer node is missing, can't collect G4 truth particles" << endl;
  }
  Bool_t tracks = kFALSE;
  Int_t nTracks = 0;
  // Iterate over tracks
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter) {
    nTracks++;
  }
  if (nTracks == 3 || nTracks == 5) {
    tracks = kTRUE;
  }

  if(tracks == kTRUE){
    return true;
  }
  else{
    return false;
  }
}


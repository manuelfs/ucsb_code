#ifndef H_EVENTHANDLER
#define H_EVENTHANDLER

#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <cfloat>
#include <stdint.h>
#include "TChain.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "pu_constants.hpp"
#include "lumi_reweighting_stand_alone.hpp"
#include "cfa.hpp"
#include "b_jet.hpp"
#include "event_number.hpp"

using namespace std;

class EventHandler : public cfA{
public:
  EventHandler(const std::string &, const bool, const double=1.0, const bool=false);
  void Skim(const std::string &, const int=-1, const int=-1);
  void MakePlots13Tev(const std::string &, int Nentries);
  void MakePlots(const std::string &);

  void SetScaleFactor(const double);
  void SetScaleFactor(const double, const double, const int);

  enum FailureType{
    kGood = 0u,
    kMETSig150 = 1u<<0,
    kMETSig100 = 1u<<1,
    kMETSig50 = 1u<<2,
    kMETSig30 = 1u<<3,
    kDeltaR = 1u<<4,
    kHiggsAvgMass = 1u<<5,
    kHiggsMassDiff = 1u<<6,
    k4thBTag = 1u<<7,
    k3rdBTag = 1u<<8,
    k2ndBTag = 1u<<9,
    kIsoTrackVeto = 1u<<10,
    kLeptonVeto = 1u<<11,
    kMinDeltaPhi = 1u<<12,
    kNumJets = 1u<<13,
    kTrigger = 1u<<14,
    kJSON = 1u<<15,
    kMETCleaning = 1u<<16,
    kPV = 1u<<17,
    kTChiMassCut = 1u<<18
  };

private:
  mutable std::pair<std::pair<TLorentzVector, TLorentzVector>, std::pair<TLorentzVector, TLorentzVector> > higgsBJetPairing;//caching for efficiency
  mutable std::vector<BJet> sortedBJetCache;//caching for efficiency
  mutable bool higgsPairingUpToDate, bJetsUpToDate;//cached value correct
  mutable bool betaUpToDate;
  static const double CSVTCut, CSVMCut, CSVLCut;
  double scaleFactor;
  mutable std::vector<double> beta;

  ///////////////////////////////////////////////////////////////////////////

  int GetTrueElectron(int iel, double &deltaR);
  int GetTrueMuon(int imu, double &deltaR);

  bool isVetoElectron(const unsigned int, const double pf_els_rel_iso_cut= 0.15) const;
  bool isVetoMuon(const unsigned int,  const double pf_mus_rel_iso_cut= 0.20) const;
  bool isVetoTau(const unsigned int) const;

  bool isRA4Electron(const unsigned int, const double pf_els_rel_iso_cut= 0.07) const;
  bool isRA4Muon(const unsigned int, const double pf_mus_rel_iso_cut= 0.10) const;

  int GetNumVetoLeptons() const;
  int GetNumVetoElectrons() const;
  int GetNumVetoMuons() const;
  int GetNumVetoTaus() const;
  int GetVetoElectron(int nth_highest_pt=1) const;
  int GetVetoMuon(int nth_highest_pt=1) const;
  float GetVetoLeptonPt(int nth_highest_pt, bool & isElectron) const;
  float GetVetoLeptonPhi(int nth_highest_pt=1) const;
  float GetVetoLeptonDeltaPhi(int nth_highest_pt=1) const;
  float GetVetoLeptonMt(int nth_highest_pt=1) const;

  int GetNumRA4Leptons() const;
  int GetNumRA4Electrons(const double pf_mus_rel_iso_cut= 0.20) const;
  int GetNumRA4Muons(const double pf_mus_rel_iso_cut= 0.10) const;
  int GetRA4Electron(int nth_highest_pt=1, const double pf_mus_rel_iso_cut= 0.07) const;
  int GetRA4Muon(int nth_highest_pt=1, const double pf_mus_rel_iso_cut= 0.10) const;
  float GetRA4LeptonPt(int nth_highest_pt, bool & isElectron) const;
  float GetRA4LeptonPhi(int nth_highest_pt=1) const;
  float GetRA4LeptonDeltaPhi(int nth_highest_pt=1) const;
  float GetRA4LeptonMt(int nth_highest_pt=1) const;

  ///////////////////////////////////////////////////////////////////////////


  int GetcfAVersion() const;

  void GetEntry(const unsigned int);

  void GetBeta(const std::string which="beta") const;

  double GetPUWeight(reweight::LumiReWeighting &) const;

  bool PassesPVCut() const;
  bool PassesMETCleaningCut() const;
  bool PassesTriggerCut() const;
  bool PassesNumJetsCut() const;
  bool Passes2CSVTCut() const;
  bool PassesJet2PtCut() const;
  bool PassesMinDeltaPhiCut() const;
  bool PassesMETSig30Cut() const;
  bool PassesLeptonVetoCut() const;
  bool PassesIsoTrackVetoCut() const;
  bool PassesBTaggingCut() const;
  bool PassesHiggsAvgMassCut() const;
  bool PassesHiggsMassDiffCut() const;
  bool PassesHiggsMassCut() const;
  bool PassesDRCut() const;
  bool PassesInvertedDRCut() const;
  bool PassesMETSig50Cut() const;
  bool PassesMETSig80Cut() const;
  bool PassesMETSig100Cut() const;
  bool PassesMETSig150Cut() const;

  bool PassesSingleLeptonCut() const;
  bool PassesJSONCut() const;

  uint_least32_t GetCutFailCode() const;

  bool PassesRegionACut() const;
  bool PassesRegionBCut() const;
  bool PassesRegionC3bCut() const;
  bool PassesRegionD3bCut() const;
  bool PassesRegionC2bCut() const;
  bool PassesRegionD2bCut() const;

  bool PassesInvertedDRRegionACut() const;
  bool PassesInvertedDRRegionBCut() const;
  bool PassesInvertedDRRegionC3bCut() const;
  bool PassesInvertedDRRegionD3bCut() const;
  bool PassesInvertedDRRegionC2bCut() const;
  bool PassesInvertedDRRegionD2bCut() const;

  bool PassesSingleLeptonRegionACut() const;
  bool PassesSingleLeptonRegionBCut() const;
  bool PassesSingleLeptonRegionC3bCut() const;
  bool PassesSingleLeptonRegionD3bCut() const;
  bool PassesSingleLeptonRegionC2bCut() const;
  bool PassesSingleLeptonRegionD2bCut() const;

  bool PassesBadJetFilter() const;

  bool PassesTChiMassCut(int=-1, int=1) const;
  bool PassesSignalMassCut(int mGlu=-1, int mLSP=1) const;
  void GetHiggsBJetPairing() const;
  void GetSortedBJets() const;

  std::pair<double, double> GetHiggsMasses() const;
  double GetHiggsDeltaR() const;

  unsigned int GetNumLowPtPfCands(const double=20.0) const;
  double GetMETOfLowPtPfCands(const double=20.0) const;

  int GetPBNR() const;
  double GetMinDeltaPhiMET(const unsigned int) const;

  int GetNumGoodJets(double ptThresh=40) const;
  int GetNumCSVTJets() const;
  int GetNumCSVMJets() const;
  int GetNumCSVLJets() const;

  bool isGoodJet(const unsigned int, const bool=true, const double=40.0, const double=2.4, const bool=true) const;
  bool isProblemJet(const unsigned int) const;
  bool jetPassLooseID(const unsigned int) const;

  bool isIsoTrack(const unsigned int, const double=10.0) const;
  bool isQualityTrack(const unsigned int) const;

  int GetNumIsoTracks(const double=10.0) const;
  int NewGetNumIsoTracks(const double=10.0) const;

  double GetElectronRelIso(const unsigned int) const;
  double GetMuonRelIso(const unsigned int) const;

  double GetSbinWeight() const;
  double GetTopPtWeight() const;

  double GetMaxDR() const;
  double GetHT(const bool=true, const bool=false) const;
  double GetHighestCSV(const unsigned int=1) const;

  bool HasGluonSplitting() const;

  std::vector<std::pair<int,int> > GetBOrigins() const;
};

#endif

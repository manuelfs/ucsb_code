#ifndef H_CUTFLOW
#define H_CUTFLOW

#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
#include <math.h>
#include "TFile.h"
#include "TChain.h"
#include "TBranch.h"

using namespace std;

class Cutflow{
public:
  Cutflow(const vector<TFile*>);
  Cutflow(const string, const bool=false);
  void PrepareVectors();
  void LoadValues();
  void Print(const bool=false) const;
  //void PrintLatex(const bool=false) const;
  void PrintCSV(const bool=false) const;
  unsigned int numCutsTotal_;
  vector<int> unscaled_;
  vector<double> scaled_;
  vector<double> error_sq_;
  vector<double> error_;
  vector<string> cutNames_;
private:
  TChain fChain_;
  
  // Declaration of leaf types
  UInt_t          startCount_;
  UInt_t          CSVTCount_;
  UInt_t          PVCount_;
  UInt_t          TriggerCount_;
  UInt_t          numJetsCount_;
  UInt_t          jet2PtCount_;
  UInt_t          minDeltaPhiCount_;
  UInt_t          METSig30Count_;
  UInt_t          METCleaningCount_;
  UInt_t          leptonVetoCount_;
  UInt_t          isoTrackVetoCount_;
  UInt_t          bTagCount_;
  UInt_t          higgsCount_;
  UInt_t          DRCount_;
  UInt_t          METSig50Count_;
  UInt_t          METSig100Count_;
  UInt_t          METSig150Count_;
  Double_t        startCountWeighted_;
  Double_t        PVCountWeighted_;
  Double_t        METCleaningCountWeighted_;
  Double_t        TriggerCountWeighted_;
  Double_t        numJetsCountWeighted_;
  Double_t        CSVTCountWeighted_;
  Double_t        jet2PtCountWeighted_;
  Double_t        minDeltaPhiCountWeighted_;
  Double_t        METSig30CountWeighted_;
  Double_t        leptonVetoCountWeighted_;
  Double_t        isoTrackVetoCountWeighted_;
  Double_t        bTagCountWeighted_;
  Double_t        higgsCountWeighted_;
  Double_t        DRCountWeighted_;
  Double_t        METSig50CountWeighted_;
  Double_t        METSig100CountWeighted_;
  Double_t        METSig150CountWeighted_;
  // List of branches
  TBranch        *b_startCount_;   //!
  TBranch        *b_CSVTCount_;   //!
  TBranch        *b_PVCount_;   //!
  TBranch        *b_TriggerCount_;   //!
  TBranch        *b_numJetsCount_;   //!
  TBranch        *b_jet2PtCount_;   //!
  TBranch        *b_minDeltaPhiCount_;   //!
  TBranch        *b_METSig30Count_;   //!
  TBranch        *b_METCleaningCount_;   //!
  TBranch        *b_leptonVetoCount_;   //!
  TBranch        *b_isoTrackVetoCount_;   //!
  TBranch        *b_bTagCount_;   //!
  TBranch        *b_higgsCount_;   //!
  TBranch        *b_DRCount_;   //!
  TBranch        *b_METSig50Count_;   //!
  TBranch        *b_METSig100Count_;   //!
  TBranch        *b_METSig150Count_;   //!
  TBranch        *b_startCountWeighted_;   //!
  TBranch        *b_PVCountWeighted_;   //!
  TBranch        *b_METCleaningCountWeighted_;   //!
  TBranch        *b_TriggerCountWeighted_;   //!
  TBranch        *b_numJetsCountWeighted_;   //!
  TBranch        *b_CSVTCountWeighted_;   //!
  TBranch        *b_jet2PtCountWeighted_;   //!
  TBranch        *b_minDeltaPhiCountWeighted_;   //!
  TBranch        *b_METSig30CountWeighted_;   //!
  TBranch        *b_leptonVetoCountWeighted_;   //!
  TBranch        *b_isoTrackVetoCountWeighted_;   //!
  TBranch        *b_bTagCountWeighted_;   //!
  TBranch        *b_higgsCountWeighted_;   //!
  TBranch        *b_DRCountWeighted_;   //!
  TBranch        *b_METSig50CountWeighted_;   //!
  TBranch        *b_METSig100CountWeighted_;   //!
  TBranch        *b_METSig150CountWeighted_;   //!
};

#endif

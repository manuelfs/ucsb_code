/* Loops over the cutflow chains in each output file, sums the event counts 
   and weights after each cut, stores them in vectors*/
#include "cutflow.hpp"

using namespace std;

Cutflow::Cutflow(const vector<TFile*> files) : 
  numCutsTotal_(13),
  unscaled_(numCutsTotal_,0),
  scaled_(numCutsTotal_,0.0),
  error_sq_(numCutsTotal_,0.0),
  error_(numCutsTotal_,0.0),
  fChain_("cutFlow")

{
  for (unsigned int f = 0; f < files.size(); f++) fChain_.Add(files[f]->GetName());
  PrepareVectors();
  LoadValues();
}

Cutflow::Cutflow(const string in_filename, bool filelist)  : 
  numCutsTotal_(13),
  unscaled_(numCutsTotal_,0),
  scaled_(numCutsTotal_,0.0),
  error_sq_(numCutsTotal_,0.0),
  error_(numCutsTotal_,0.0),
  fChain_("cutFlow")
{
  char path[1024]="";
  if(filelist){
    FILE *fp(fopen(in_filename.c_str(), "r"));
    while(!feof(fp)){
      char line[1024]="";
      int scanner(fscanf(fp, "%s", line));
      sprintf(path, "%s", line);
      fChain_.Add(path);
      if(false && scanner==0){
      }
    }
    fclose(fp);
  }
  else{
    sprintf(path, "%s", in_filename.c_str());
    fChain_.Add(path);
  } 
        
  PrepareVectors();
  LoadValues();
}


void Cutflow::PrepareVectors() {
  /*char *cutNames[numCutsTotal_]  = {"Start", // 0 (cutflow index)
    "Primary Vertex", // 1
    "Jet 2 pT > 50 GeV", // 2
    "METsig>30", // 3
    "MET Cleaning", // 4
    "Trigger", // 5
    "nJets = 4 || 5", // 6
    "Min. Delta Phi >0.3", // 7
    "Lepton Veto", // 8
    "Isolated Track Veto", // 9
    "4 b-Jets", // 10
    "Higgs mass window cuts", // 11
    "Max. delta R<2.2"}; // 12*/
  cutNames_.push_back("Start"); // 0 (cutflow index)
  cutNames_.push_back("PV"); // 1
  cutNames_.push_back("$\\pt_{j_2}>50$ GeV"); // 2
  cutNames_.push_back("Cleaning"); // 3
  cutNames_.push_back("$\\metsig>30$"); // 4
  cutNames_.push_back("Presel."); // 5 (trigger)
  cutNames_.push_back("nJets $= 4 || 5$"); // 6
  cutNames_.push_back("$\\mdp>0.3$"); // 7
  cutNames_.push_back("Lepton Veto"); // 8
  cutNames_.push_back("IsoTk. Veto"); // 9
  cutNames_.push_back("4 b-Jets"); // 10
  cutNames_.push_back("$\\left<m_{bb}\\right>,\\Delta m_{bb}$"); // 11
  cutNames_.push_back("$\\mdR<2.2$"); // 12

}
void Cutflow::LoadValues() {
  fChain_.SetBranchAddress("startCount", &startCount_, &b_startCount_);
  fChain_.SetBranchAddress("CSVTCount", &CSVTCount_, &b_CSVTCount_);
  fChain_.SetBranchAddress("PVCount", &PVCount_, &b_PVCount_);
  fChain_.SetBranchAddress("TriggerCount", &TriggerCount_, &b_TriggerCount_);
  fChain_.SetBranchAddress("numJetsCount", &numJetsCount_, &b_numJetsCount_);
  fChain_.SetBranchAddress("jet2PtCount", &jet2PtCount_, &b_jet2PtCount_);
  fChain_.SetBranchAddress("minDeltaPhiCount", &minDeltaPhiCount_, &b_minDeltaPhiCount_);
  fChain_.SetBranchAddress("METSig30Count", &METSig30Count_, &b_METSig30Count_);
  fChain_.SetBranchAddress("METCleaningCount", &METCleaningCount_, &b_METCleaningCount_);
  fChain_.SetBranchAddress("leptonVetoCount", &leptonVetoCount_, &b_leptonVetoCount_);
  fChain_.SetBranchAddress("isoTrackVetoCount", &isoTrackVetoCount_, &b_isoTrackVetoCount_);
  fChain_.SetBranchAddress("bTagCount", &bTagCount_, &b_bTagCount_);
  fChain_.SetBranchAddress("higgsCount", &higgsCount_, &b_higgsCount_);
  fChain_.SetBranchAddress("DRCount", &DRCount_, &b_DRCount_);
  fChain_.SetBranchAddress("METSig50Count", &METSig50Count_, &b_METSig50Count_);
  fChain_.SetBranchAddress("METSig100Count", &METSig100Count_, &b_METSig100Count_);
  fChain_.SetBranchAddress("METSig150Count", &METSig150Count_, &b_METSig150Count_);
  fChain_.SetBranchAddress("startCountWeighted", &startCountWeighted_, &b_startCountWeighted_);
  fChain_.SetBranchAddress("PVCountWeighted", &PVCountWeighted_, &b_PVCountWeighted_);
  fChain_.SetBranchAddress("METCleaningCountWeighted", &METCleaningCountWeighted_, &b_METCleaningCountWeighted_);
  fChain_.SetBranchAddress("TriggerCountWeighted", &TriggerCountWeighted_, &b_TriggerCountWeighted_);
  fChain_.SetBranchAddress("numJetsCountWeighted", &numJetsCountWeighted_, &b_numJetsCountWeighted_);
  fChain_.SetBranchAddress("CSVTCountWeighted", &CSVTCountWeighted_, &b_CSVTCountWeighted_);
  fChain_.SetBranchAddress("jet2PtCountWeighted", &jet2PtCountWeighted_, &b_jet2PtCountWeighted_);
  fChain_.SetBranchAddress("minDeltaPhiCountWeighted", &minDeltaPhiCountWeighted_, &b_minDeltaPhiCountWeighted_);
  fChain_.SetBranchAddress("METSig30CountWeighted", &METSig30CountWeighted_, &b_METSig30CountWeighted_);
  fChain_.SetBranchAddress("leptonVetoCountWeighted", &leptonVetoCountWeighted_, &b_leptonVetoCountWeighted_);
  fChain_.SetBranchAddress("isoTrackVetoCountWeighted", &isoTrackVetoCountWeighted_, &b_isoTrackVetoCountWeighted_);
  fChain_.SetBranchAddress("bTagCountWeighted", &bTagCountWeighted_, &b_bTagCountWeighted_);
  fChain_.SetBranchAddress("higgsCountWeighted", &higgsCountWeighted_, &b_higgsCountWeighted_);
  fChain_.SetBranchAddress("DRCountWeighted", &DRCountWeighted_, &b_DRCountWeighted_);
  fChain_.SetBranchAddress("METSig50CountWeighted", &METSig50CountWeighted_, &b_METSig50CountWeighted_);
  fChain_.SetBranchAddress("METSig100CountWeighted", &METSig100CountWeighted_, &b_METSig100CountWeighted_);
  fChain_.SetBranchAddress("METSig150CountWeighted", &METSig150CountWeighted_, &b_METSig150CountWeighted_);
                
  for (unsigned int entry = 0; entry < fChain_.GetEntries(); entry++) {
    fChain_.GetEntry(entry);
    unscaled_[0]+=startCount_;
    unscaled_[1]+=PVCount_;
    unscaled_[2]+=jet2PtCount_;
    unscaled_[3]+=METSig30Count_;
    unscaled_[4]+=METCleaningCount_;
    unscaled_[5]+=TriggerCount_;
    unscaled_[6]+=numJetsCount_;
    unscaled_[7]+=minDeltaPhiCount_;
    unscaled_[8]+=leptonVetoCount_;
    unscaled_[9]+=isoTrackVetoCount_;
    unscaled_[10]+=bTagCount_;
    unscaled_[11]+=higgsCount_;
    unscaled_[12]+=DRCount_;
    scaled_[0]+=startCountWeighted_;
    scaled_[1]+=PVCountWeighted_;
    scaled_[2]+=jet2PtCountWeighted_;
    scaled_[3]+=METSig30CountWeighted_;
    scaled_[4]+=METCleaningCountWeighted_;
    scaled_[5]+=TriggerCountWeighted_;
    scaled_[6]+=numJetsCountWeighted_;
    scaled_[7]+=minDeltaPhiCountWeighted_;
    scaled_[8]+=leptonVetoCountWeighted_;
    scaled_[9]+=isoTrackVetoCountWeighted_;
    scaled_[10]+=bTagCountWeighted_;
    scaled_[11]+=higgsCountWeighted_;
    scaled_[12]+=DRCountWeighted_;
    for (unsigned int e = 0; e < error_sq_.size(); e++) {
      if (unscaled_[e]>0) error_sq_[e]+=(scaled_[e]*scaled_[e])/(double)unscaled_[e];
    }
  }
  for (unsigned int e = 0; e < error_sq_.size(); e++) {
    error_[e]=sqrt(error_sq_[e]);
  }
} 

void Cutflow::Print(const bool unweighted) const{
  for (unsigned int c = 0; c<numCutsTotal_; c++) {
    if (unweighted) printf("%s: %d\n",cutNames_[c].c_str(),unscaled_[c]);
    else printf("%s: %f +/- %f\n",cutNames_[c].c_str(),scaled_[c],error_[c]);
  }
}

void Cutflow::PrintCSV(const bool unweighted) const{ // comma-separated values (for Excel)
  for (unsigned int c = 0; c<numCutsTotal_; c++) {
    if (unweighted) printf("%s, %d\n",cutNames_[c].c_str(),unscaled_[c]);
    else printf("%s, %f, %f\n",cutNames_[c].c_str(),scaled_[c],error_[c]);
  }
}


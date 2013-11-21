#include "TString.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
//#include "macros/PlotUtils.cc"
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>

#define Ntrees 1

using namespace std;
using std::cout;
using std::endl;


 


int GetcfAVersion(string sampleName){
  size_t pos(sampleName.rfind("_v"));
  if(pos!=std::string::npos && pos<(sampleName.size()-2)){
    std::istringstream iss(sampleName.substr(pos+2));
    int version(0);
    iss >> version;
    return version;
  }else{
    return 0;
  }
}

bool jetPassLooseID(const unsigned int ijet, std::vector<float> *jets_AK5PF_energy,
		    std::vector<float> *jets_AK5PF_corrFactorRaw, std::vector<float> *jets_AK5PFclean_mu_Mult,
		    std::vector<float> *jets_AK5PF_neutral_Mult, std::vector<float> *jets_AK5PF_chg_Mult, 
		    std::vector<float> *jets_AK5PF_neutralHadE, std::vector<float>* jets_AK5PF_chgHadE, 
		    std::vector<float> *jets_AK5PFclean_neutralEmE, std::vector<float> *jets_AK5PF_chgEmE,
		    std::vector<float> *jets_AK5PF_eta){

  const double jetenergy = jets_AK5PF_energy->at(ijet) * jets_AK5PF_corrFactorRaw->at(ijet);
  const int numConst = static_cast<int>(jets_AK5PFclean_mu_Mult->at(ijet)+jets_AK5PF_neutral_Mult->at(ijet)+
					jets_AK5PF_chg_Mult->at(ijet));
  
  if (jetenergy>0.0) {
    if (jets_AK5PF_neutralHadE->at(ijet) /jetenergy <= 0.99
	&& jets_AK5PFclean_neutralEmE->at(ijet) / jetenergy <= 0.99
	&& numConst >= 2
	&& ( fabs(jets_AK5PF_eta->at(ijet))>=2.4
	     || (fabs(jets_AK5PF_eta->at(ijet))<2.4 && jets_AK5PF_chgHadE->at(ijet)/jetenergy>0))
	&& ( fabs(jets_AK5PF_eta->at(ijet))>=2.4
	     || (fabs(jets_AK5PF_eta->at(ijet))<2.4 && jets_AK5PF_chgEmE->at(ijet)/jetenergy<0.99))
	&& ( fabs(jets_AK5PF_eta->at(ijet))>=2.4
	     || (fabs(jets_AK5PF_eta->at(ijet))<2.4 && jets_AK5PF_chg_Mult->at(ijet)>0))){
      return true;
    }
  }
  return false;
}

bool isGoodJet(const unsigned int ijet, const bool jetid, const double ptThresh, const double etaThresh, 
	       std::vector<float> *jets_AK5PF_energy, std::vector<float> *jets_AK5PF_corrFactorRaw, 
	       std::vector<float> *jets_AK5PFclean_mu_Mult, std::vector<float> *jets_AK5PF_neutral_Mult,
	       std::vector<float> *jets_AK5PF_chg_Mult, std::vector<float> *jets_AK5PF_neutralHadE,
	       std::vector<float>* jets_AK5PF_chgHadE, 
	       std::vector<float> *jets_AK5PFclean_neutralEmE, std::vector<float> *jets_AK5PF_chgEmE,
	       std::vector<float>*jets_AK5PF_pt,
	       std::vector<float> *jets_AK5PF_eta) {
  if(jets_AK5PF_pt->at(ijet)<ptThresh || fabs(jets_AK5PF_eta->at(ijet))>etaThresh) return false;
  if(jetid && !jetPassLooseID(ijet, jets_AK5PF_energy, jets_AK5PF_corrFactorRaw, jets_AK5PFclean_mu_Mult,
			      jets_AK5PF_neutral_Mult, jets_AK5PF_chg_Mult,
			      jets_AK5PF_neutralHadE, jets_AK5PF_chgHadE,
			      jets_AK5PFclean_neutralEmE, jets_AK5PF_chgEmE,
			      jets_AK5PF_eta)) return false;
      
     return true;
}

void compare13tev(int iVariable = 0){
  // Variables
  int Nbins[] = {50,50,50};
  double limits[][2] = {{0,100},{-6,6},{0,15}};
  TString Variables[] = {"Njets_AK5PF","jets_AK5PF_eta", "NGoodJets"};
   
 

  // Ntuples
  TString folder("/net/cms2/cms2r0/cfA/"), name, hname[Ntrees];
  TString treeNames[] = {"TTbar_TuneZ2star_13TeV-pythia6-tauola_Summer13dr53X-PU45bx25_START53_V19D-v2_AODSIM_UCSB1933_v71"};
  TString tagNames[] = {"ttbar"};
  TChain tree[Ntrees]; 
  for(int itree(0); itree < Ntrees; itree++){
    name = folder; name += treeNames[itree]; name += "/*.root/configurableAnalysis/eventB";
    tree[itree].Add(name);
  }
  

  TCanvas can("can","8 TeV Vs 13/14 TeV comparison");
  TH1F* htree[Ntrees];
  for(int itree(0); itree < Ntrees; itree++){
    hname[itree] = tagNames[itree]; hname[itree] += itree;
    htree[itree] = new TH1F(hname[itree],Variables[iVariable],Nbins[iVariable], limits[iVariable][0], limits[iVariable][1]);
    htree[itree]->SetLineWidth(2);
    //tree[itree].Project(hname[itree], Variables[iVariable]);
 
    //Branches
    tree[itree].SetBranchStatus("*",0);
    tree[itree].SetBranchStatus("Njets_AK5PF",1);
    tree[itree].SetBranchStatus("jets_AK5PF_pt",1);
    tree[itree].SetBranchStatus("jets_AK5PF_eta",1);

    UInt_t Njets_AK5PF;
    tree[itree].SetBranchAddress("Njets_AK5PF", &Njets_AK5PF);
    std::vector<float> *jets_AK5PF_pt(0);
    tree[itree].SetBranchAddress("jets_AK5PF_pt", &jets_AK5PF_pt);
    std::vector<float> *jets_AK5PF_eta(0);
    tree[itree].SetBranchAddress("jets_AK5PF_eta", &jets_AK5PF_eta);
    std::vector<float> *jets_AK5PF_energy(0);
    tree[itree].SetBranchAddress("jets_AK5PF_energy", &jets_AK5PF_energy);
    std::vector<float> *jets_AK5PF_corrFactorRaw(0);
    tree[itree].SetBranchAddress("jets_AK5PF_corrFactorRaw", &jets_AK5PF_corrFactorRaw);
    std::vector<float> *jets_AK5PF_neutralHadE(0);
    tree[itree].SetBranchAddress("jets_AK5PF_neutralHadE", &jets_AK5PF_neutralHadE);
    std::vector<float> *jets_AK5PFclean_neutralEmE(0);
    tree[itree].SetBranchAddress("jets_AK5PFclean_neutralEmE", &jets_AK5PFclean_neutralEmE);
    std::vector<float> *jets_AK5PFclean_mu_Mult(0);  
    tree[itree].SetBranchAddress("jets_AK5PFclean_mu_Mult", &jets_AK5PFclean_mu_Mult);
    std::vector<float> *jets_AK5PF_neutral_Mult(0);
    tree[itree].SetBranchAddress("jets_AK5PF_neutral_Mult", &jets_AK5PF_neutral_Mult);
    std::vector<float> *jets_AK5PF_chg_Mult(0);
    tree[itree].SetBranchAddress("jets_AK5PF_chg_Mult", &jets_AK5PF_chg_Mult);
    std::vector<float> *jets_AK5PF_chgHadE(0);
    tree[itree].SetBranchAddress("jets_AK5PF_chgHadE", &jets_AK5PF_chgHadE);
    std::vector<float> *jets_AK5PF_chgEmE(0);
    tree[itree].SetBranchAddress("jets_AK5PF_chgEmE", &jets_AK5PF_chgEmE);



    int Nentries = tree[itree].GetEntries();
    Nentries = 10000;
    for(int entry(0); entry < Nentries; entry++){
      tree[itree].GetEntry(entry);
      int NGoodJets=0;
      for(unsigned int ijet(0); ijet < jets_AK5PF_pt->size(); ijet++ ){
      //      htree[itree]->Fill(jets_AK5PF_pt->size());    
	if (isGoodJet(ijet, true, 40, 2.4, jets_AK5PF_energy, jets_AK5PF_corrFactorRaw, 
		      jets_AK5PFclean_mu_Mult, jets_AK5PF_neutral_Mult, jets_AK5PF_chg_Mult, 
		      jets_AK5PF_neutralHadE, jets_AK5PF_chgHadE, jets_AK5PFclean_neutralEmE, 
		      jets_AK5PF_chgEmE, jets_AK5PF_pt, jets_AK5PF_eta))
	  NGoodJets++;
      } // Loop over all jets
      htree[itree]->Fill(NGoodJets);  
    } // Loop over all events
    htree[itree]->Draw();
  
}
  
  TString pName = "plots/"; pName += Variables[iVariable]; pName += ".pdf";
  can.SaveAs(pName);
  for(int itree(0); itree < Ntrees; itree++) htree[itree]->Delete();



}


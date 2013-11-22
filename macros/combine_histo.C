#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TTree.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TKey.h"
#include "TList.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TLatex.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>

#define NFiles 6

using namespace std;
using std::cout;
using std::endl;


void combine_histo(){
 
  //Files
  TString FileNames[] = {"raw_plots_and_values/SMS-T1tttt_2J_mGo-845to3000_mLSP-1to1355_TuneZ2star_14TeV-madgraph-tauola_Summer12-START53_V7C_FSIM_PU_S12-v1_AODSIM_UCSB1949reshuf_v71.root",
			 "raw_plots_and_values/SMS-MadGraph_Pythia6Zstar_8TeV_T1tttt_2J_mGo-1100to1400_mLSP-525to1000_25GeVX25GeV_Binning_Summer12-START52_V9_FSIM-v2_AODSIM_UCSB1739reshuf_v68.root",
			 "raw_plots_and_values/TTbar_TuneZ2star_13TeV-pythia6-tauola_Summer13dr53X-PU45bx25_START53_V19D-v2_AODSIM_UCSB1933_v71.root",
			 "raw_plots_and_values/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71.root",
			 "raw_plots_and_values/TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71.root",
			 "raw_plots_and_values/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71.root"};
  
  TString tagNames[] = {"T1tttt_14TeV", "T1tttt_8TeV", "ttbar_13TeV", "ttbar_ll_8TeV", "ttbar_hh_8TeV", "ttbar_lh_8TeV"}, Pname;
  TFile *file[NFiles];
  for(int iFiles(0); iFiles < NFiles; iFiles++)
    file[iFiles] = new TFile(FileNames[iFiles]);

  int colors[] = {2, 4, 3, 1, 6, 9};
  TH1F hFile[NFiles];
  TCanvas can;
  //Loop over all variables  
  for(int obj(0); obj < file[0]->GetListOfKeys()->GetSize(); ++obj){
    const std::string obj_name(file[0]->GetListOfKeys()->At(obj)->GetName());

    for(int iFiles(0); iFiles < NFiles; iFiles++){
	hFile[iFiles] = *(static_cast<TH1F*>(file[iFiles]->GetKey(obj_name.c_str(),1)->ReadObj()));
	hFile[iFiles].SetLineColor(colors[iFiles]);
	if(iFiles==0) hFile[iFiles].Draw("h");
	else hFile[iFiles].Draw("same h");
      } //Loop over all files
    Pname = "plots/"; Pname += obj_name; Pname += ".pdf";
    can.SaveAs(Pname);
  } //Loop over all variables

  for(int iFiles(0); iFiles < NFiles; iFiles++){
    file[iFiles]->Close();
    file[iFiles]->Delete();
  }
}
    

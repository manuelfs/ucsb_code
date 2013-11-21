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

#define NFiles 2

using namespace std;
using std::cout;
using std::endl;


void combine_histo(){
 
  //Files
  TString FileNames[] = {"raw_plots_and_values/TTbar_TuneZ2star_13TeV-pythia6-tauola_Summer13dr53X-PU45bx25_START53_V19D-v2_AODSIM_UCSB1933_v71.root", 
			 "raw_plots_and_values/SMS-T1tttt_2J_mGo-825_mLSP-225_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V19-v1_AODSIM_UCSB1788reshuf_v68.root",
			 "raw_plots_and_values/SMS-T1tttt_2J_mGo-845to3000_mLSP-1to1355_TuneZ2star_14TeV-madgraph-tauola_Summer12-START53_V7C_FSIM_PU_S12-v1_AODSIM_UCSB1949reshuf_v71.root"};
  TString tagNames[] = {"ttbar", "T1tttt_mGo_825", "T1tttt_mGo_845_3000"};
  TFile *file[NFiles];
  for(int iFiles(0); iFiles < NFiles; iFiles++)
    file[iFiles] = new TFile(FileNames[iFiles]);

  int colors[] = {2, 4, 3};
  TH1F hFile[NFiles];
  TCanvas can;
  //Loop over all variables  
  //for(int obj(0); obj < file[0]->GetListOfKeys()->GetSize(); ++obj){
  for(int obj(0); obj < 1; ++obj){
    cout<<"Getting histo name "<<obj<<endl;
    const std::string obj_name(file[0]->GetListOfKeys()->At(obj)->GetName());
    cout<<"Plotting histogram "<<obj_name<<endl;

    for(int iFiles(0); iFiles < NFiles; iFiles++){
	hFile[iFiles] = *(static_cast<TH1F*>(file[iFiles]->GetKey(obj_name.c_str(),1)->ReadObj()));
	hFile[iFiles].SetLineColor(colors[iFiles]);
	if(iFiles==0) hFile[iFiles].Draw("h");
	else hFile[iFiles].Draw("same h");
      } //Loop over all variables 

  } //Loop over all files
  can.SaveAs("test.pdf");

  for(int iFiles(0); iFiles < NFiles; iFiles++){
    file[iFiles]->Close();
    file[iFiles]->Delete();
  }
}
    

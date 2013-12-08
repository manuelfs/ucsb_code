#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TTree.h"
#include "TColor.h"
#include "TStyle.h"
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

TString RoundNumber(double n, int e, double d=1);


void combine_histo(){
 
  //Files
  TString FileNames[] = { "raw_plots_and_values/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71.root",
			  "raw_plots_and_values/TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71.root",
			  "raw_plots_and_values/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71.root",
			  "raw_plots_and_values/TTbar_TuneZ2star_13TeV-pythia6-tauola_Summer13dr53X-PU45bx25_START53_V19D-v2_AODSIM_UCSB1933_v71.root",
			  "raw_plots_and_values/SMS-MadGraph_Pythia6Zstar_8TeV_T1tttt_2J_mGo-1100to1400_mLSP-525to1000_25GeVX25GeV_Binning_Summer12-START52_V9_FSIM-v2_AODSIM_UCSB1739reshuf_v68_1150_800_.root",
			  "raw_plots_and_values/SMS-T1tttt_2J_mGo-845to3000_mLSP-1to1355_TuneZ2star_14TeV-madgraph-tauola_Summer12-START53_V7C_FSIM_PU_S12-v1_AODSIM_UCSB1949reshuf_v71_1145_800.root"};
  
  TString tagNames[] = {"ttbar_ll_8TeV", "ttbar_hh_8TeV", "ttbar_lh_8TeV", "ttbar_13TeV","T1tttt_8TeV","T1tttt_14TeV"}, Pname;
  TString legNames[] = {"ttbar_ll_8TeV", "ttbar_hh_8TeV", "t#bar{t} @ 8 TeV ", "t#bar{t} @ 13 TeV ","T1tttt(1150,800) @ 8 TeV ","T1tttt(1145,800) @ 14 TeV "};
  //TString xTitles[] = {"Number of good jets", "H_{T} (GeV)", "E_{T,miss} (GeV)", "LOpt", "NLOpt", "NNLOpt" };
  // TString yTitles[] = {"Entries", "Entries/(100 GeV)", "Entries/(50 GeV)", "Entries", "Entries", "Entries"};
  //TString texNames[] = {"$\\left<N_j \\right>$", "$\\left<H_T \\right>$", "$\\left<E_{T, {\\rm miss}} \\right>$"}; 
  TFile *file[NFiles];
  for(int iFiles(0); iFiles < NFiles; iFiles++)
    file[iFiles] = new TFile(FileNames[iFiles]);
  
  int colors[] = {kSpring-9, kSpring-9, kSpring-9, kGreen+2, kRed-9, kRed+1};
  int styles[] = {3345, 3345, 3345, 0, 3354, 0};

  ofstream texFile; texFile.open("txtrgbt/Averages.tex");

  gStyle->SetHatchesLineWidth(2);
  gStyle->SetOptStat(0);
  TH1F hFile[NFiles];
  TCanvas can;
  //Loop over all variables  
  TString xTitle = "", VarName, texName;
  TString yTitle = "";
  TString Title = "";
  for(int obj(0); obj < file[0]->GetListOfKeys()->GetSize(); ++obj){
    const std::string obj_name(file[0]->GetListOfKeys()->At(obj)->GetName());
    VarName = obj_name;
    if(VarName.Contains("NumGoodJets")) {xTitle = "Number of Good Jets "; texName = "$\\left<N_j \\right>$";}
    if(VarName.Contains("HT")) xTitle = "H_{T} ";
    if(VarName.Contains("MET")) xTitle = "E_{T,miss} ";
    if(VarName.Contains("NumLeptonVeto")) xTitle = "Number of Veto Leptons";
    if(VarName.Contains("NumLeptonRA4")) xTitle = "Number of RA4 Leptons";
    if(VarName.Contains("MT")) xTitle = "M_{T} ";
    if(VarName.Contains("pTRA41")) xTitle = "p_{T} of First RA4 Lepton ";
    if(VarName.Contains("pTVeto1")) xTitle = "p_{T} of First Veto Lepton ";
    if(VarName.Contains("pTRA42")) xTitle = "p_{T} of Second RA4 Lepton ";
    if(VarName.Contains("pTVeto2")) xTitle = "p_{T} of Second Veto Lepton ";
    if(VarName.Contains("pT_")) Title = "p_{T} ";
    if(VarName.Contains("HT")||VarName.Contains("MET")||VarName.Contains("MT")||VarName.Contains("pTVeto")||VarName.
       Contains("pTRA4")||VarName.Contains("pT_")) xTitle+="(GeV)";				   
    //if(VarName.Contains("_1l"))  Title+= "Single Lepton ";
    // if(VarName.Contains("3jets_"))  Title+= " 3 Jets ";
    // if(VarName.Contains("l_3jets"))  Title+= " and 3 Jets ";
    // if(VarName.Contains("l_4jets"))  Title+= " and 4 jets ";
    // if(VarName.Contains("l_5jets"))  Title+= " and 5 jets ";

    double maxhisto(0), means[NFiles];  
    TLegend leg(0.5,0.7,0.89,0.89);
    leg.SetTextSize(0.04); leg.SetFillColor(0); leg.SetBorderSize(0);
    leg.SetTextFont(132);
    for(int iFiles(0); iFiles < NFiles; iFiles++){
      hFile[iFiles] = *(static_cast<TH1F*>(file[iFiles]->GetKey(obj_name.c_str(),1)->ReadObj()));
      hFile[iFiles].SetLineColor(colors[iFiles]);
      hFile[iFiles].SetLineWidth(2);
      hFile[iFiles].SetFillStyle(styles[iFiles]);
      hFile[iFiles].SetFillColor(colors[iFiles]);
      if(iFiles==2){ 
	hFile[iFiles].Add(&hFile[0],.25);
	hFile[iFiles].Add(&hFile[1],1 );
      }
      Title=hFile[iFiles].GetTitle();
      if(VarName.Contains("jets")){ 
	Title+=", p_{T}^{thresh} = ";
	TString VarName2 = VarName;
	VarName2.Remove(0, VarName2.Sizeof()-2);
	Title+=VarName2;
	Title+= " GeV";
      }    
      yTitle = "Entries ";
      if(VarName.Contains("HT")||VarName.Contains("MET")||VarName.Contains("MT")||VarName.Contains("pTVeto")||VarName.Contains("pTRA4")||VarName.Contains("pT_")){				   
	yTitle+="/ (";
	yTitle+= RoundNumber(hFile[iFiles].GetBinWidth(1), 0);
	yTitle+=" GeV)";
	  }
      hFile[iFiles].SetTitle(Title);
      hFile[iFiles].SetXTitle(xTitle);
      hFile[iFiles].SetYTitle(yTitle);
      hFile[iFiles].Scale(1000./hFile[iFiles].Integral());    
      cout<<"The Mean of "<< tagNames[iFiles]<<" is "<<hFile[iFiles].GetMean()<<endl;
      if(hFile[iFiles].GetMaximum() > maxhisto && iFiles>=2) maxhisto = hFile[iFiles].GetMaximum();
      if(iFiles>=2) leg.AddEntry(&hFile[iFiles], legNames[iFiles]);
      means[iFiles] = hFile[iFiles].GetMean();
    } //Loop over all files
  
    texFile << "\\begin{tabular}{c | ccc}\n \\hline\\hline\n"<<texName<<" & 8 TeV & 13/14 TeV & $\\Delta$ (\\%) \\\\"<<endl;
    texFile << "\\hline\n $t\\bar{t}$ & "<< RoundNumber(means[2],1) << " & " << RoundNumber(means[3],1) << " & "<<
      RoundNumber((means[3]-means[2])*100, 1, means[2]) << "\\\\" << endl;
    texFile << "T1tttt & "<< RoundNumber(means[4],1) << " & " << RoundNumber(means[5],1) << " & "<<
      RoundNumber((means[5]-means[4])*100, 1, means[4]) << "\\\\ \\hline\\hline" << endl;
    texFile << "\\end{tabular}"<<endl<<endl;
    for(int iFiles(0); iFiles < NFiles; iFiles++){
      if(iFiles==2){
	hFile[iFiles].SetMaximum(1.15*maxhisto); 
	hFile[iFiles].Draw("h");
      }else hFile[iFiles].Draw("same h");
    } //Loop over all files
    hFile[3].Draw("same h");    
    leg.Draw();
    Pname = "plots/"; Pname += obj_name; Pname += ".pdf";
    can.SetLogy(0);    
    can.SaveAs(Pname);
    hFile[2].SetMaximum(10*maxhisto);    
    can.SetLogy(1);
    Pname.ReplaceAll(".pdf", "_log.pdf");
    can.SaveAs(Pname);
  } //Loop over all variables

  texFile.close();
  for(int iFiles(0); iFiles < NFiles; iFiles++){
    file[iFiles]->Close();
    file[iFiles]->Delete();
  }
}
    
TString RoundNumber(double n, int e, double d){
  if(d==0) return " - ";
  double neg = 1; if(n*d<0) neg = -1;
  n /= neg*d; n += 0.5*pow(10.,-e);
  int n_int = (int)n;
  int n_dec = (int)((1+n-n_int)*pow(10.,e));
  TString s_dec = ""; s_dec += n_dec; s_dec.Remove(0,1);
  TString result=""; 
  if(neg<0) result+="-";
  result+= n_int;
  if(e>0) {
    result+="."; result+=s_dec;
  }
  
  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<e-afterdot.Length(); i++)
    result += "0";
  return result;
}

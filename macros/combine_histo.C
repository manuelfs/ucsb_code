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

#define NFiles 7

using namespace std;
using std::cout;
using std::endl;

TString RoundNumber(double n, int e, double d=1);


void combine_histo(){
 
  //Files
  TString FileNames[] = { "raw_plots_and_values/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71_.root",
			  "raw_plots_and_values/TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71_.root",
			  "raw_plots_and_values/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71_.root",
			  "raw_plots_and_values/TTbar_TuneZ2star_13TeV-pythia6-tauola_Summer13dr53X-PU45bx25_START53_V19D-v2_AODSIM_UCSB1933_v71_.root",
			  "raw_plots_and_values/SMS-MadGraph_Pythia6Zstar_8TeV_T1tttt_2J_mGo-1100to1400_mLSP-525to1000_25GeVX25GeV_Binning_Summer12-START52_V9_FSIM-v2_AODSIM_UCSB1739reshuf_v68_1150_800_.root",
			  "raw_plots_and_values/SMS-T1tttt_2J_mGo-845to3000_mLSP-1to1355_TuneZ2star_14TeV-madgraph-tauola_Summer12-START53_V7C_FSIM_PU_S12-v1_AODSIM_UCSB1949reshuf_v71_1145_800.root", 
			  "raw_plots_and_values/SMS-T1tttt_2J_mGo-845to3000_mLSP-1to1355_TuneZ2star_14TeV-madgraph-tauola_Summer12-START53_V7C_FSIM_PU_S12-v1_AODSIM_UCSB1949reshuf_v71_1145_500.root"};
  
  // TString tagNames[] = {"ttbar_ll_8TeV", "ttbar_hh_8TeV", "ttbar_lh_8TeV", "ttbar_13TeV",
  // 			"T1tttt_8TeV", "T1tttt_14TeV_1145_800", "T1tttt_14TeV_1145_500"};
  TString legNames[] = {"ttbar_ll_8TeV", "ttbar_hh_8TeV", "t#bar{t} @ 8 TeV ", "t#bar{t} @ 13 TeV ",
			"T1tttt(1150,800) @ 8 TeV ","T1tttt(1145,800) @ 14 TeV ", "T1tttt(1145,500) @ 14 TeV "}, Pname;
  TFile *file[NFiles];
  for(int iFiles(0); iFiles < NFiles; iFiles++){
    //    cout<<"Opening "<<FileNames[iFiles]<<endl;
    file[iFiles] = new TFile(FileNames[iFiles]);
  }
  
  int colors[] = {kSpring-9, kSpring-9, kSpring-9, kGreen+2, kRed-9, kRed+1, kBlue+1};
  int styles[] = {3345, 3345, 3345, 0, 3354, 0, 0};

  ofstream texFile; texFile.open("txt/Averages.tex");

  gStyle->SetCanvasDefW(1000);
  gStyle->SetCanvasDefH(600);

  gStyle->SetHatchesLineWidth(2);
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.038);
  gStyle->SetPadTopMargin(0.035);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.1);
  TH1F hFile[NFiles];
  TCanvas can;
  //Loop over all variables  
  TString xTitle = "", VarName, texName;
  TString yTitle = "";
  TString Title = "";
  for(int obj(0); obj < file[0]->GetListOfKeys()->GetSize(); ++obj){
    const std::string obj_name(file[0]->GetListOfKeys()->At(obj)->GetName());
    cout<< "Doing variable "<<obj_name<<endl;
    VarName = obj_name;
    xTitle = "";
    texName = "";
    if(VarName.Contains("HT")) {xTitle = "H_{T} "; texName = "$\\bf{\\left<H_T \\right>} \\mathbf{(GeV)}$";};
    if(VarName.Contains("MET")) {xTitle = "E_{T,miss} ";  texName = "$\\bf{\\left<E_{T,miss} \\right>} \\mathbf{(GeV)}$";};
    if(VarName.Contains("MT"))  {xTitle = "M_{T} "; texName = "$\\bf{\\left<M_T \\right>} \\mathbf{(GeV)}$";};
    if(VarName.Contains("pTRA41"))  {xTitle = "p_T of First RA4 Lepton "; texName = "$\\bf{\\left<p_T \\right>} \\mathbf{(GeV)}$";};
    if(VarName.Contains("pTVeto1")) { xTitle = "p_{T} of First Veto Lepton "; texName = "$\\bf{\\left<p_T \\right>} \\mathbf{(GeV)}$";};
    if(VarName.Contains("pTRA42")) {xTitle = "p_{T} of Second RA4 Lepton "; texName = "$\\bf{\\left<p_T \\right>} \\mathbf{(GeV)}$";};
    if(VarName.Contains("pTVeto2")){ xTitle = "p_{T} of Second Veto Lepton "; texName = "$\\bf{\\left<p_T \\right>} \\mathbf{(GeV)}$";};
    if(VarName.Contains("pT_")) {xTitle = "p_{T} "; texName = "$\\bf{\\left<p_T \\right>} \\mathbf{(GeV)}$";};
    xTitle+="(GeV)";
				   
    if(VarName.Contains("NumLeptonsVeto")) {xTitle = "Number of Veto Leptons "; texName = "$\\bf{\\left<N_veto \\right>} $";};
    if(VarName.Contains("NumLeptonsRA4")) {xTitle = "Number of RA4 Leptons"; texName ="$\\bf{\\left<N_RA4 \\right>} $"; };
    if(VarName.Contains("NumGoodJets")) {xTitle = "Number of Good Jets "; texName = "$\\bf{\\left<N_j \\right>} $";};
    if(VarName.Contains("Iso")) {xTitle = "Relative isolation"; texName = "$\\bf{\\left<Iso_{\\rm rel} \\right>} $";};
    if(VarName.Contains("PU")) {xTitle = "Number of true interactions"; texName = "$\\bf{\\left<N_{\\rm Vert} \\right>} $";};

    double maxhisto(0), means[NFiles];  
    TLegend leg(0.45,0.65,0.96,0.90);
    leg.SetTextSize(0.06); leg.SetFillColor(0); leg.SetBorderSize(0);
    leg.SetTextFont(132);
    for(int iFiles(0); iFiles < NFiles; iFiles++){
      //cout<< "Getting keys from file "<< iFiles<<", pointer "<<file[iFiles]<<endl;
      hFile[iFiles] = *(static_cast<TH1F*>(file[iFiles]->GetKey(obj_name.c_str(),1)->ReadObj()));
      // cout<< hFile[iFiles].GetTitle()<<endl;
      hFile[iFiles].SetLineColor(colors[iFiles]);
      hFile[iFiles].SetLineWidth(2);
      hFile[iFiles].SetFillStyle(styles[iFiles]);
      hFile[iFiles].SetFillColor(colors[iFiles]);
      if(iFiles==2){ 
	hFile[iFiles].Add(&hFile[0],.25);
	hFile[iFiles].Add(&hFile[1],1 );
      }
      Title=hFile[iFiles].GetTitle();
      if(VarName.Contains("ets")){ 
	Title+=", p_{T}^{thresh} = ";
	TString VarName2 = VarName;
	VarName2.Remove(0, VarName2.Sizeof()-3);
	Title+=VarName2;
	Title+= " GeV";
      }    
      yTitle = "Entries ";
      if(VarName.Contains("HT")||VarName.Contains("MET")||VarName.Contains("MT")||VarName.Contains("pTVeto")||VarName.Contains("pTRA4")||VarName.Contains("pT_")){				   
	yTitle+="/ (";
	yTitle+= RoundNumber(hFile[iFiles].GetBinWidth(1), 0);
	yTitle+=" GeV)";
      }
      if(VarName.Contains("Iso")){				   
	yTitle+="/ (";
	yTitle+= RoundNumber(hFile[iFiles].GetBinWidth(1), 2);
	yTitle+=")";
      }
      hFile[iFiles].SetTitle(Title);
      hFile[iFiles].SetXTitle(xTitle);
      hFile[iFiles].SetYTitle(yTitle);
      //hFile[iFiles].SetTextSize(0.06);            // Set global text size
      hFile[iFiles].SetTitleSize(0.05,"xy");     // Set the 2 axes title size
      hFile[iFiles].SetLabelSize(0.05,"xy");     // Set the 2 axes label size
      if(iFiles>=2) hFile[iFiles].Scale(1000./hFile[iFiles].Integral());    
      //cout<<"The Mean of "<< tagNames[iFiles]<<" is "<<hFile[iFiles].GetMean()<<endl;
      if(hFile[iFiles].GetMaximum() > maxhisto && iFiles>=2) maxhisto = hFile[iFiles].GetMaximum();
      if(iFiles>=2) leg.AddEntry(&hFile[iFiles], legNames[iFiles]);
      means[iFiles] = hFile[iFiles].GetMean();
    } //Loop over all files
    //cout<<"Writing table"<<endl;   
    int digits = 1;
    if(VarName.Contains("Iso")) digits = 3;
    texFile << obj_name << endl;
    texFile << "\\begin{tabular}{c | ccc}\n \\hline\\hline\n"<<texName<<" & 8 TeV & 13/14 TeV & $\\Delta$ (\\%) \\\\"<<endl;
    texFile << "\\hline\n $t\\bar{t}$ & "<< RoundNumber(means[2],digits) << " & " << RoundNumber(means[3],digits) << " & "<<
      RoundNumber((means[3]-means[2])*100, 1, means[2]) << "\\\\" << endl;
    texFile << "T1tttt & "<< RoundNumber(means[4],digits) << " & " << RoundNumber(means[5],digits) << " & "<<
      RoundNumber((means[5]-means[4])*100, 1, means[4]) << "\\\\ \\hline\\hline" << endl;
    texFile << "\\end{tabular}"<<endl<<endl;
    for(int iFiles(0); iFiles < NFiles; iFiles++){
      if(iFiles==2){
	hFile[iFiles].SetMaximum(1.2*maxhisto); 
	if(VarName.Contains("PU")) hFile[iFiles].SetMaximum(100); 
	hFile[iFiles].Draw("h");
      }else hFile[iFiles].Draw("same h");
    } //Loop over all files
    hFile[3].Draw("same h");    
    leg.Draw();
    Pname = "plots/"; Pname += obj_name; Pname += ".pdf";
    can.SetLogy(0);    
    can.SaveAs(Pname);
    hFile[2].SetMaximum(20*maxhisto);    
    if(VarName.Contains("PU")) hFile[2].SetMaximum(100); 
    can.SetLogy(1);
    hFile[2].SetMinimum(0.05);
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

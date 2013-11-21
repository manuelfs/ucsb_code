#include <vector>
#include <string>
#include <cmath>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TKey.h"
#include "style.hpp"

TH1D GetTH1D(std::vector<TFile*> &file, const std::string &obj_name){
  if(file.size()){
    TH1D h= *(static_cast<TH1D*>(file.at(0)->GetKey(obj_name.c_str(),1)->ReadObj()));
    for(unsigned int i(1); i<file.size(); ++i){
      h=h+ *(static_cast<TH1D*>(file.at(i)->GetKey(obj_name.c_str(),1)->ReadObj()));
    }
    return h;
  }else{
    TH1D h;
    return h;
  }
}

TH1D GetSOverSqrtB(TH1D sig, TH1D back){
  if(sig.GetNbinsX()!=back.GetNbinsX()) return sig;
  sig.GetYaxis()->SetTitle("S/#sqrt{B}");
  for(int i(1); i<=sig.GetNbinsX(); ++i){
    if(back.GetBinContent(i)>0.0 && sig.GetBinContent(i)>0.0){
      sig.SetBinContent(i,sig.GetBinContent(i)/sqrt(back.GetBinContent(i)));
    }else{
      sig.SetBinContent(i,0.0);
    }
  }
  return sig;
}

void ConvertToIntegral(TH1D &h){
  const int N(h.GetNbinsX());
  for(int i(0); i<=N+1; ++i){
    h.SetBinContent(i,h.Integral(i,N+1));
  }
}

int main(){
  SetStyle();
  std::vector<std::string> plotName(0);
  plotName.push_back("sorb_MET");
  plotName.push_back("sorb_METSig");
  plotName.push_back("sorb_METSig2012");
  plotName.push_back("sorb_METOverSqrtHT");

  for(unsigned int i(0); i<plotName.size(); ++i){
    std::vector<TFile*> SM(0), signal400(0), signal250(0);
    
    SM.push_back(new TFile("raw_plots_and_values/QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM_UCSB1513_v66_Adam.root","read"));
    SM.push_back(new TFile("raw_plots_and_values/QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1603_v66_Adam.root","read"));
    SM.push_back(new TFile("raw_plots_and_values/QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_v3_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1609_v66_Adam.root","read"));
    SM.push_back(new TFile("raw_plots_and_values/QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1659_v67_Adam.root","read"));
    SM.push_back(new TFile("raw_plots_and_values/QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1663_v67_Adam.root","read"));
    SM.push_back(new TFile("raw_plots_and_values/QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1559_v66_Adam.root","read"));
    SM.push_back(new TFile("raw_plots_and_values/QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1664_v67_Adam.root","read"));
    SM.push_back(new TFile("raw_plots_and_values/QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1665_v67_Adam.root","read"));
    SM.push_back(new TFile("raw_plots_and_values/QCD_Pt-1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1585_v66_Adam.root","read"));

    SM.push_back(new TFile("raw_plots_and_values/TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66_Jack.root","read"));

    SM.push_back(new TFile("raw_plots_and_values/TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1606_v66_Jack.root","read"));

    SM.push_back(new TFile("raw_plots_and_values/TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1605_v66_Adam.root","read"));
    SM.push_back(new TFile("raw_plots_and_values/TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1604_v66_Adam.root","read"));

    SM.push_back(new TFile("raw_plots_and_values/WbbJetsToLNu_Massive_TuneZ2star_8TeV-madgraph-pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1677_v67_Adam.root","read"));
  
    SM.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66_Adam.root","read"));
    SM.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1607_v66_Adam.root","read"));
    SM.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1524_v66_Adam.root","read"));
    SM.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1594_v66_Adam.root","read"));
    SM.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1523_v66_Adam.root","read"));
    SM.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1602_v66_Adam.root","read"));

    SM.push_back(new TFile("raw_plots_and_values/TTH_Inclusive_M-125_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM+_UCSB1783_v68_Jack.root","read"));

    signal400.push_back(new TFile("raw_plots_and_values/configurableAnalysis_TChihh_400_1_Adam.root","read"));

    signal250.push_back(new TFile("raw_plots_and_values/configurableAnalysis_TChihh_250_1_Adam.root","read"));
    TCanvas canvas;

    TH1D h_SM(GetTH1D(SM, plotName.at(i)));
    TH1D h_signal250(GetTH1D(signal250, plotName.at(i)));
    TH1D h_signal400(GetTH1D(signal400, plotName.at(i)));

    h_SM.SetStats(0);
    h_signal250.SetStats(0);
    h_signal400.SetStats(0);

    h_SM.Draw("hist");
    canvas.Print(("plots/SOverSqrtB_h_SM_"+plotName.at(i)+".pdf").c_str());
    h_signal250.Draw("hist");
    canvas.Print(("plots/SOverSqrtB_h_signal250_"+plotName.at(i)+".pdf").c_str());
    h_signal400.Draw("hist");
    canvas.Print(("plots/SOverSqrtB_h_signal400_"+plotName.at(i)+".pdf").c_str());

    ConvertToIntegral(h_SM);
    ConvertToIntegral(h_signal250);
    ConvertToIntegral(h_signal400);

    h_SM.Draw("hist");
    canvas.Print(("plots/SOverSqrtB_h_SM_Integral_"+plotName.at(i)+".pdf").c_str());
    h_signal250.Draw("hist");
    canvas.Print(("plots/SOverSqrtB_h_signal250_Integral_"+plotName.at(i)+".pdf").c_str());
    h_signal400.Draw("hist");
    canvas.Print(("plots/SOverSqrtB_h_signal400_Integral_"+plotName.at(i)+".pdf").c_str());

    TH1D sorb250(GetSOverSqrtB(h_signal250, h_SM));
    TH1D sorb400(GetSOverSqrtB(h_signal400, h_SM));

    sorb250.SetLineColor(2);
    sorb400.SetLineColor(1);

    double max(sorb250.GetMaximum());
    if(sorb400.GetMaximum()>max){
      max=sorb400.GetMaximum();
    }
    max*=1.1;
    sorb250.SetMaximum(max);
    sorb400.SetMaximum(max);

    sorb250.Draw("hist");
    sorb400.Draw("histsame");

    canvas.Print(("plots/SOverSqrtB_"+plotName.at(i)+".pdf").c_str());
  }
}

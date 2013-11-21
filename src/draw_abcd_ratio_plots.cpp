#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>
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
#include "style.hpp"

/*const int c_QCD(TColor::GetColor(127,255,255));
  const int c_TTJets_FullLept(TColor::GetColor(63,63,255));
  const int c_TTJets_SemiLept(TColor::GetColor(127,127,255));
  const int c_TTJets_Hadronic(TColor::GetColor(191,191,255));
  const int c_TTV(TColor::GetColor(127,255,127));
  const int c_Wbb(TColor::GetColor(255,127,255));
  const int c_ZJets(TColor::GetColor(255,127,127));
  const int c_TTH(TColor::GetColor(127,127,127));
  const int c_realData(TColor::GetColor(0,0,0));
  const int c_backgroundMC(TColor::GetColor(0,0,255));
  const int c_signal400(TColor::GetColor(255,0,0));
  const int c_signal250(TColor::GetColor(0,255,0));*/

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

TH2D GetTH2D(std::vector<TFile*> &file, const std::string &obj_name){
  if(file.size()){
    TH2D h= *(static_cast<TH2D*>(file.at(0)->GetKey(obj_name.c_str(),1)->ReadObj()));
    for(unsigned int i(1); i<file.size(); ++i){
      h=h+ *(static_cast<TH2D*>(file.at(i)->GetKey(obj_name.c_str(),1)->ReadObj()));
    }
    return h;
  }else{
    TH2D h;
    return h;
  }
}

TH2D operator+(const TH2D &h1, const TH2D &h2){
  if(h1.GetNbinsX()!=h2.GetNbinsX() || h1.GetNbinsY()!=h2.GetNbinsY()){
    TH2D h;
    return h;
  }else{
    TH2D h(h1);
    for(int x(1); x<=h1.GetNbinsX(); ++x){
      for(int y(1); y<h1.GetNbinsY(); ++y){
        h.SetBinContent(x,y,h1.GetBinContent(x,y)+h2.GetBinContent(x,y));
      }
    }
    return h;
  }
}

int main(int argc, char *argv[]){
  TH1::SetDefaultSumw2(true);
  SetStyle();
  if(argc>1){
    std::vector<TFile*> QCD(0), TTJets_FullLept(0), TTJets_SemiLept(0), TTJets_Hadronic(0), TTV(0), Wbb(0), ZJets(0), TTH(0), signal400(0), signal250(0), realData(0);
    
    QCD.push_back(new TFile("raw_plots_and_values/QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM_UCSB1654_v67.root","read"));
    QCD.push_back(new TFile("raw_plots_and_values/QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1657_v67.root","read"));
    QCD.push_back(new TFile("raw_plots_and_values/QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_v3_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1658_v67.root","read"));
    QCD.push_back(new TFile("raw_plots_and_values/QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1659_v67.root","read"));
    QCD.push_back(new TFile("raw_plots_and_values/QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1663_v67.root","read"));
    QCD.push_back(new TFile("raw_plots_and_values/QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1660_v67.root","read"));
    QCD.push_back(new TFile("raw_plots_and_values/QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1664_v67.root","read"));
    QCD.push_back(new TFile("raw_plots_and_values/QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1665_v67.root","read"));
    QCD.push_back(new TFile("raw_plots_and_values/QCD_Pt-1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1666_v67.root","read"));

    TTJets_FullLept.push_back(new TFile("raw_plots_and_values/TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66.root","read"));

    TTJets_SemiLept.push_back(new TFile("raw_plots_and_values/TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1606_v66.root","read"));

    TTJets_Hadronic.push_back(new TFile("raw_plots_and_values/TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1613_v66.root","read"));

    TTV.push_back(new TFile("raw_plots_and_values/TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1605_v66.root","read"));
    TTV.push_back(new TFile("raw_plots_and_values/TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1604_v66.root","read"));

    Wbb.push_back(new TFile("raw_plots_and_values/WbbJetsToLNu_Massive_TuneZ2star_8TeV-madgraph-pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1677_v67.root","read"));
  
    ZJets.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1607_v66.root","read"));
    ZJets.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66.root","read"));
    ZJets.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1524_v66.root","read"));
    ZJets.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1594_v66.root","read"));
    ZJets.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1523_v66_Adam.root","read"));
    ZJets.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1602_v66_Adam.root","read"));

    TTH.push_back(new TFile("raw_plots_and_values/TTH_Inclusive_M-125_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM+_UCSB1783_v68.root","read"));

    realData.push_back(new TFile("raw_plots_and_values/Run2012A.root","read"));
    realData.push_back(new TFile("raw_plots_and_values/Run2012B.root","read"));
    realData.push_back(new TFile("raw_plots_and_values/Run2012C.root","read"));
    realData.push_back(new TFile("raw_plots_and_values/Run2012D.root","read"));

    signal400.push_back(new TFile("raw_plots_and_values/SMS-HbbHbb_mHiggsino-400_mLSP-1_8TeV-Pythia6Z_jgsmith-SMS-HbbHbb_mHiggsino-400_mLSP-1_8TeV-Pythia6Z-26439e701cfb9736f297615863e915f9_USER_UCSB1812_v69.root","read"));

    signal250.push_back(new TFile("raw_plots_and_values/SMS-HbbHbb_mHiggsino-250_mLSP-1_8TeV-Pythia6Z_jgsmith-SMS-HbbHbb_mHiggsino-250_mLSP-1_8TeV-Pythia6Z-26439e701cfb9736f297615863e915f9_USER_UCSB1808_v69.root","read"));

    //TFile outFile((argv[1]+std::string(".root")).c_str(), "recreate");

    TH1D A_QCD(GetTH1D(QCD,"xx_metSig_A"));
    TH1D A_TTJets_FullLept(GetTH1D(TTJets_FullLept,"xx_metSig_A"));
    TH1D A_TTJets_SemiLept(GetTH1D(TTJets_SemiLept,"xx_metSig_A"));
    TH1D A_TTJets_Hadronic(GetTH1D(TTJets_Hadronic,"xx_metSig_A"));
    TH1D A_TTV(GetTH1D(TTV,"xx_metSig_A"));
    TH1D A_Wbb(GetTH1D(Wbb,"xx_metSig_A"));
    TH1D A_ZJets(GetTH1D(ZJets,"xx_metSig_A"));
    TH1D A_TTH(GetTH1D(TTH,"xx_metSig_A"));
    TH1D A_signal400(GetTH1D(signal400,"xx_metSig_A"));
    TH1D A_signal250(GetTH1D(signal250,"xx_metSig_A"));
    TH1D A_SM(A_QCD+A_TTJets_FullLept+A_TTJets_SemiLept+A_TTJets_Hadronic+A_TTV+A_Wbb+A_ZJets+A_TTH);

    TH1D B_QCD(GetTH1D(QCD,"xx_metSig_B"));
    TH1D B_TTJets_FullLept(GetTH1D(TTJets_FullLept,"xx_metSig_B"));
    TH1D B_TTJets_SemiLept(GetTH1D(TTJets_SemiLept,"xx_metSig_B"));
    TH1D B_TTJets_Hadronic(GetTH1D(TTJets_Hadronic,"xx_metSig_B"));
    TH1D B_TTV(GetTH1D(TTV,"xx_metSig_B"));
    TH1D B_Wbb(GetTH1D(Wbb,"xx_metSig_B"));
    TH1D B_ZJets(GetTH1D(ZJets,"xx_metSig_B"));
    TH1D B_TTH(GetTH1D(TTH,"xx_metSig_B"));
    TH1D B_signal400(GetTH1D(signal400,"xx_metSig_B"));
    TH1D B_signal250(GetTH1D(signal250,"xx_metSig_B"));
    TH1D B_SM(B_QCD+B_TTJets_FullLept+B_TTJets_SemiLept+B_TTJets_Hadronic+B_TTV+B_Wbb+B_ZJets+B_TTH);

    TH1D C3b_QCD(GetTH1D(QCD,"xx_metSig_C3b"));
    TH1D C3b_TTJets_FullLept(GetTH1D(TTJets_FullLept,"xx_metSig_C3b"));
    TH1D C3b_TTJets_SemiLept(GetTH1D(TTJets_SemiLept,"xx_metSig_C3b"));
    TH1D C3b_TTJets_Hadronic(GetTH1D(TTJets_Hadronic,"xx_metSig_C3b"));
    TH1D C3b_TTV(GetTH1D(TTV,"xx_metSig_C3b"));
    TH1D C3b_Wbb(GetTH1D(Wbb,"xx_metSig_C3b"));
    TH1D C3b_ZJets(GetTH1D(ZJets,"xx_metSig_C3b"));
    TH1D C3b_TTH(GetTH1D(TTH,"xx_metSig_C3b"));
    TH1D C3b_signal400(GetTH1D(signal400,"xx_metSig_C3b"));
    TH1D C3b_signal250(GetTH1D(signal250,"xx_metSig_C3b"));
    TH1D C3b_SM(C3b_QCD+C3b_TTJets_FullLept+C3b_TTJets_SemiLept+C3b_TTJets_Hadronic+C3b_TTV+C3b_Wbb+C3b_ZJets+C3b_TTH);

    TH1D D3b_QCD(GetTH1D(QCD,"xx_metSig_D3b"));
    TH1D D3b_TTJets_FullLept(GetTH1D(TTJets_FullLept,"xx_metSig_D3b"));
    TH1D D3b_TTJets_SemiLept(GetTH1D(TTJets_SemiLept,"xx_metSig_D3b"));
    TH1D D3b_TTJets_Hadronic(GetTH1D(TTJets_Hadronic,"xx_metSig_D3b"));
    TH1D D3b_TTV(GetTH1D(TTV,"xx_metSig_D3b"));
    TH1D D3b_Wbb(GetTH1D(Wbb,"xx_metSig_D3b"));
    TH1D D3b_ZJets(GetTH1D(ZJets,"xx_metSig_D3b"));
    TH1D D3b_TTH(GetTH1D(TTH,"xx_metSig_D3b"));
    TH1D D3b_signal400(GetTH1D(signal400,"xx_metSig_D3b"));
    TH1D D3b_signal250(GetTH1D(signal250,"xx_metSig_D3b"));
    TH1D D3b_SM(D3b_QCD+D3b_TTJets_FullLept+D3b_TTJets_SemiLept+D3b_TTJets_Hadronic+D3b_TTV+D3b_Wbb+D3b_ZJets+D3b_TTH);

    TH1D C2b_QCD(GetTH1D(QCD,"xx_metSig_C2b"));
    TH1D C2b_TTJets_FullLept(GetTH1D(TTJets_FullLept,"xx_metSig_C2b"));
    TH1D C2b_TTJets_SemiLept(GetTH1D(TTJets_SemiLept,"xx_metSig_C2b"));
    TH1D C2b_TTJets_Hadronic(GetTH1D(TTJets_Hadronic,"xx_metSig_C2b"));
    TH1D C2b_TTV(GetTH1D(TTV,"xx_metSig_C2b"));
    TH1D C2b_Wbb(GetTH1D(Wbb,"xx_metSig_C2b"));
    TH1D C2b_ZJets(GetTH1D(ZJets,"xx_metSig_C2b"));
    TH1D C2b_TTH(GetTH1D(TTH,"xx_metSig_C2b"));
    TH1D C2b_signal400(GetTH1D(signal400,"xx_metSig_C2b"));
    TH1D C2b_signal250(GetTH1D(signal250,"xx_metSig_C2b"));
    TH1D C2b_SM(C2b_QCD+C2b_TTJets_FullLept+C2b_TTJets_SemiLept+C2b_TTJets_Hadronic+C2b_TTV+C2b_Wbb+C2b_ZJets+C2b_TTH);

    TH1D D2b_QCD(GetTH1D(QCD,"xx_metSig_D2b"));
    TH1D D2b_TTJets_FullLept(GetTH1D(TTJets_FullLept,"xx_metSig_D2b"));
    TH1D D2b_TTJets_SemiLept(GetTH1D(TTJets_SemiLept,"xx_metSig_D2b"));
    TH1D D2b_TTJets_Hadronic(GetTH1D(TTJets_Hadronic,"xx_metSig_D2b"));
    TH1D D2b_TTV(GetTH1D(TTV,"xx_metSig_D2b"));
    TH1D D2b_Wbb(GetTH1D(Wbb,"xx_metSig_D2b"));
    TH1D D2b_ZJets(GetTH1D(ZJets,"xx_metSig_D2b"));
    TH1D D2b_TTH(GetTH1D(TTH,"xx_metSig_D2b"));
    TH1D D2b_signal400(GetTH1D(signal400,"xx_metSig_D2b"));
    TH1D D2b_signal250(GetTH1D(signal250,"xx_metSig_D2b"));
    TH1D D2b_SM(D2b_QCD+D2b_TTJets_FullLept+D2b_TTJets_SemiLept+D2b_TTJets_Hadronic+D2b_TTV+D2b_Wbb+D2b_ZJets+D2b_TTH);

    TCanvas canvas;

    TH1D r4b_SM(A_SM/B_SM);
    TH1D r3b_SM(C3b_SM/D3b_SM);
    TH1D r2b_SM(C2b_SM/D2b_SM);
    r4b_SM.SetMarkerColor(1);
    r3b_SM.SetMarkerColor(2);
    r2b_SM.SetMarkerColor(3);
    r4b_SM.SetLineColor(1);
    r3b_SM.SetLineColor(2);
    r2b_SM.SetLineColor(3);
    r4b_SM.SetMarkerStyle(20);
    r3b_SM.SetMarkerStyle(20);
    r2b_SM.SetMarkerStyle(20);
    r4b_SM.SetStats(0);
    r3b_SM.SetStats(0);
    r2b_SM.SetStats(0);

    r2b_SM.Draw("e1p");
    r3b_SM.Draw("e1psame");
    r4b_SM.Draw("e1psame");
    
    TLegend legend_SM(0.7,0.5,0.95,0.85);
    legend_SM.SetFillColor(4000);
    legend_SM.AddEntry(&r4b_SM, "4b", "l");
    legend_SM.AddEntry(&r3b_SM, "3b", "l");
    legend_SM.AddEntry(&r2b_SM, "2b", "l");
    legend_SM.Draw("same");

    canvas.Print((argv[1]+std::string("_SM_SigOverSB")+".pdf").c_str());
    canvas.SetLogy(1);
    canvas.Print((argv[1]+std::string("_SM_SigOverSB")+"_log.pdf").c_str());
    canvas.SetLogy(0);

    TH1D r4b_signal400(A_signal400/B_signal400);
    TH1D r3b_signal400(C3b_signal400/D3b_signal400);
    TH1D r2b_signal400(C2b_signal400/D2b_signal400);
    r4b_signal400.SetMarkerColor(1);
    r3b_signal400.SetMarkerColor(2);
    r2b_signal400.SetMarkerColor(3);
    r4b_signal400.SetLineColor(1);
    r3b_signal400.SetLineColor(2);
    r2b_signal400.SetLineColor(3);
    r4b_signal400.SetMarkerStyle(20);
    r3b_signal400.SetMarkerStyle(20);
    r2b_signal400.SetMarkerStyle(20);
    r4b_signal400.SetStats(0);
    r3b_signal400.SetStats(0);
    r2b_signal400.SetStats(0);

    r2b_signal400.Draw("e1p");
    r3b_signal400.Draw("e1psame");
    r4b_signal400.Draw("e1psame");
    
    TLegend legend_signal400(0.7,0.5,0.95,0.85);
    legend_signal400.SetFillColor(4000);
    legend_signal400.AddEntry(&r4b_signal400, "4b", "l");
    legend_signal400.AddEntry(&r3b_signal400, "3b", "l");
    legend_signal400.AddEntry(&r2b_signal400, "2b", "l");
    legend_signal400.Draw("same");

    canvas.Print((argv[1]+std::string("_signal400_SigOverSB")+".pdf").c_str());
    canvas.SetLogy(1);
    canvas.Print((argv[1]+std::string("_signal400_SigOverSB")+"_log.pdf").c_str());
    canvas.SetLogy(0);

    TH1D r4b_signal250(A_signal250/B_signal250);
    TH1D r3b_signal250(C3b_signal250/D3b_signal250);
    TH1D r2b_signal250(C2b_signal250/D2b_signal250);
    r4b_signal250.SetMarkerColor(1);
    r3b_signal250.SetMarkerColor(2);
    r2b_signal250.SetMarkerColor(3);
    r4b_signal250.SetLineColor(1);
    r3b_signal250.SetLineColor(2);
    r2b_signal250.SetLineColor(3);
    r4b_signal250.SetMarkerStyle(20);
    r3b_signal250.SetMarkerStyle(20);
    r2b_signal250.SetMarkerStyle(20);
    r4b_signal250.SetStats(0);
    r3b_signal250.SetStats(0);
    r2b_signal250.SetStats(0);

    r2b_signal250.Draw("e1p");
    r3b_signal250.Draw("e1psame");
    r4b_signal250.Draw("e1psame");
    
    TLegend legend_signal250(0.7,0.5,0.95,0.85);
    legend_signal250.SetFillColor(4000);
    legend_signal250.AddEntry(&r4b_signal250, "4b", "l");
    legend_signal250.AddEntry(&r3b_signal250, "3b", "l");
    legend_signal250.AddEntry(&r2b_signal250, "2b", "l");
    legend_signal250.Draw("same");

    canvas.Print((argv[1]+std::string("_signal250_SigOverSB")+".pdf").c_str());
    canvas.SetLogy(1);
    canvas.Print((argv[1]+std::string("_signal250_SigOverSB")+"_log.pdf").c_str());
    canvas.SetLogy(0);
    //outFile.Close();
  }
}

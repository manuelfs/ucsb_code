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

const int c_QCD(kOrange-4);
const int c_TTJets_FullLept(kBlue-7);
const int c_TTJets_SemiLept(kBlue);
const int c_TTJets_Hadronic(kBlue+2);
const int c_TTH(kGreen+1);
const int c_TTV(kGreen-4);
const int c_VH(kMagenta+1);
const int c_VV(kMagenta-4);
const int c_V(kCyan);
const int c_singleT(kGray);
const int c_realData(kBlack);
const int c_backgroundMC(kSpring+6);
const int c_signal400(kRed+1);
const int c_signal250(kRed-4);

std::string GetScatArg(const TH2D &h, std::string s=""){
  const double theMax(h.GetBinContent(h.GetMaximumBin()));
  const double scaling(theMax>0.0?500.0/theMax:1.0);
  s+="scat=";
  char scaleNum[256]="";
  sprintf(scaleNum,"%f",scaling);
  s+=scaleNum;
  return s;
}

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
    std::vector<TFile*> QCD(0), TTJets_FullLept(0), TTJets_SemiLept(0),
      TTJets_Hadronic(0), TTV(0), V(0), TTH(0), VV(0), VH(0), singleT(0),
      signal400(0), signal250(0), realData(0);
    
    realData.push_back(new TFile("raw_plots_and_values/MET_Run2012A-13Jul2012-v1_AOD_UCSB1852_v71_SyncSkim.root","read"));
    realData.push_back(new TFile("raw_plots_and_values/MET_Run2012B-13Jul2012-v1_AOD_UCSB1853_v71_SyncSkim.root","read"));
    realData.push_back(new TFile("raw_plots_and_values/MET_Run2012C-24Aug2012-v1_AOD_UCSB1854_v71_SyncSkim.root","read"));
    realData.push_back(new TFile("raw_plots_and_values/MET_Run2012C-PromptReco-v2_AOD_UCSB1867_v71_SyncSkim.root","read"));
    realData.push_back(new TFile("raw_plots_and_values/MET_Run2012D-PromptReco-v1_AOD_UCSB1869_v71_SyncSkim.root","read"));
    realData.push_back(new TFile("raw_plots_and_values/MET_Run2012D-PromptReco-v1_AOD_UCSB1870_v71_SyncSkim.root","read"));
    signal250.push_back(new TFile("raw_plots_and_values/SMS-TChiHH_2b2b_2J_mChargino-250_mLSP-1_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872_v71_SyncSkim.root","read"));
    signal400.push_back(new TFile("raw_plots_and_values/SMS-TChiHH_2b2b_2J_mChargino-400_mLSP-1_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1871_v71_SyncSkim.root","read"));
    QCD.push_back(new TFile("raw_plots_and_values/BJets_HT-250To500_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1893_v71_SyncSkim.root","read"));
    QCD.push_back(new TFile("raw_plots_and_values/BJets_HT-500To1000_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1894_v71_SyncSkim.root","read"));
    QCD.push_back(new TFile("raw_plots_and_values/BJets_HT-1000ToInf_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1895_v71_SyncSkim.root","read"));
    singleT.push_back(new TFile("raw_plots_and_values/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1864_v71_SyncSkim.root","read"));
    singleT.push_back(new TFile("raw_plots_and_values/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1865_v71_SyncSkim.root","read"));
    singleT.push_back(new TFile("raw_plots_and_values/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1866_v71_SyncSkim.root","read"));
    singleT.push_back(new TFile("raw_plots_and_values/T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1860_v71_SyncSkim.root","read"));
    singleT.push_back(new TFile("raw_plots_and_values/T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1861_v71_SyncSkim.root","read"));
    singleT.push_back(new TFile("raw_plots_and_values/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1862_v71_SyncSkim.root","read"));
    TTH.push_back(new TFile("raw_plots_and_values/TTH_HToBB_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1855_v71_SyncSkim.root","read"));
    TTJets_FullLept.push_back(new TFile("raw_plots_and_values/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71_SyncSkim.root","read"));
    TTJets_Hadronic.push_back(new TFile("raw_plots_and_values/TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71_SyncSkim.root","read"));
    TTJets_SemiLept.push_back(new TFile("raw_plots_and_values/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71_SyncSkim.root","read"));
    TTV.push_back(new TFile("raw_plots_and_values/TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1857_v71_SyncSkim.root","read"));
    TTV.push_back(new TFile("raw_plots_and_values/TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1856_v71_SyncSkim.root","read"));
    V.push_back(new TFile("raw_plots_and_values/W2JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1877_v71_SyncSkim.root","read"));
    V.push_back(new TFile("raw_plots_and_values/W3JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1878_v71_SyncSkim.root","read"));
    V.push_back(new TFile("raw_plots_and_values/W4JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1879_v71_SyncSkim.root","read"));
    V.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1887_v71_SyncSkim.root","read"));
    V.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1889_v71_SyncSkim.root","read"));
    V.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1888_v71_SyncSkim.root","read"));
    V.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1891_v71_SyncSkim.root","read"));
    V.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1890_v71_SyncSkim.root","read"));
    VH.push_back(new TFile("raw_plots_and_values/WH_WToLNu_HToBB_M-125_8TeV-powheg-herwigpp_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1858_v71_SyncSkim.root","read"));
    //VH.push_back(new TFile("raw_plots_and_values/ZH_ZToBB_HToBB_M-125_8TeV-powheg-herwigpp_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1868_v71_SyncSkim.root","read"));
    VV.push_back(new TFile("raw_plots_and_values/WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1874_v71_SyncSkim.root","read"));
    VV.push_back(new TFile("raw_plots_and_values/ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1876_v71_SyncSkim.root","read"));
    
    TFile outFile((argv[1]+std::string(".root")).c_str(), "recreate");

    for(int obj(0); obj<QCD.at(0)->GetListOfKeys()->GetSize(); ++obj){
      const std::string obj_name(QCD.at(0)->GetListOfKeys()->At(obj)->GetName());
      const double signalScale((obj_name.find("xx_")!=std::string::npos)?1.0:20.0);

      if(QCD.at(0)->GetKey(obj_name.c_str(),1)->ReadObj()->IsA()==TH1D::Class()){
        TH1D h_QCD(GetTH1D(QCD,obj_name));
        TH1D h_TTJets_FullLept(GetTH1D(TTJets_FullLept,obj_name));
        TH1D h_TTJets_SemiLept(GetTH1D(TTJets_SemiLept,obj_name));
        TH1D h_TTJets_Hadronic(GetTH1D(TTJets_Hadronic,obj_name));
        TH1D h_TTV(GetTH1D(TTV,obj_name));
        TH1D h_TTH(GetTH1D(TTH,obj_name));
        TH1D h_VH(GetTH1D(VH,obj_name));
        TH1D h_V(GetTH1D(V,obj_name));
        TH1D h_VV(GetTH1D(VV,obj_name));
        TH1D h_singleT(GetTH1D(singleT,obj_name));
        TH1D h_realData(GetTH1D(realData,obj_name));
        TH1D h_signal400(GetTH1D(signal400,obj_name));
        TH1D h_signal250(GetTH1D(signal250,obj_name));

        h_QCD.SetFillColor(c_QCD);
        h_TTJets_FullLept.SetFillColor(c_TTJets_FullLept);
        h_TTJets_SemiLept.SetFillColor(c_TTJets_SemiLept);
        h_TTJets_Hadronic.SetFillColor(c_TTJets_Hadronic);
        h_TTV.SetFillColor(c_TTV);
        h_TTH.SetFillColor(c_TTH);
        h_VH.SetFillColor(c_VH);
        h_VV.SetFillColor(c_VV);
        h_V.SetFillColor(c_V);
        h_singleT.SetFillColor(c_singleT);
        h_realData.SetFillColor(c_realData);
        h_signal400.SetFillColor(4000);
        h_signal250.SetFillColor(4000);
        h_signal400.SetFillStyle(0);
        h_signal250.SetFillStyle(0);

        h_QCD.SetLineColor(c_QCD);
        h_TTJets_FullLept.SetLineColor(c_TTJets_FullLept);
        h_TTJets_SemiLept.SetLineColor(c_TTJets_SemiLept);
        h_TTJets_Hadronic.SetLineColor(c_TTJets_Hadronic);
        h_TTV.SetLineColor(c_TTV);
        h_TTH.SetLineColor(c_TTH);
        h_VH.SetLineColor(c_VH);
        h_VV.SetLineColor(c_VV);
        h_V.SetLineColor(c_V);
        h_singleT.SetLineColor(c_singleT);
        h_realData.SetLineColor(c_realData);
        h_signal400.SetLineColor(c_signal400);
        h_signal250.SetLineColor(c_signal250);
        h_signal400.SetLineWidth(3);
        h_signal250.SetLineWidth(3);

        const double scaling=1.0;//h_realData.Integral("width")/(h_QCD.Integral("width")+h_TTJets_FullLept.Integral("width")+h_TTJets_SemiLept.Integral("width")+h_TTJets_Hadronic.Integral("width")+h_TTV.Integral("width")+h_Wbb.Integral("width")+h_ZJets.Integral("width")+h_TTH.Integral("width"));
        h_QCD.Scale(scaling);
        h_TTJets_FullLept.Scale(scaling);
        h_TTJets_SemiLept.Scale(scaling);
        h_TTJets_Hadronic.Scale(scaling);
        h_TTV.Scale(scaling);
        h_TTH.Scale(scaling);
        h_VH.Scale(scaling);
        h_VV.Scale(scaling);
        h_V.Scale(scaling);
        h_singleT.Scale(scaling);
        h_signal400.Scale(scaling);
        h_signal250.Scale(scaling);
        
        h_signal400.Scale(signalScale);
        h_signal250.Scale(signalScale);
          
        std::string title((std::string(h_QCD.GetTitle())+";"
                           +h_QCD.GetXaxis()->GetTitle())+";"
                          +h_QCD.GetYaxis()->GetTitle());
        THStack stack(obj_name.c_str(),title.c_str());
        stack.Add(&h_TTJets_Hadronic);
        stack.Add(&h_TTJets_SemiLept);
        stack.Add(&h_TTJets_FullLept);
        stack.Add(&h_TTH);
        stack.Add(&h_TTV);
        stack.Add(&h_VH);
        stack.Add(&h_VV);
        stack.Add(&h_V);
        stack.Add(&h_singleT);
        stack.Add(&h_QCD);

        int dataMaxBin=h_realData.GetMaximumBin();
        double themax=h_realData.GetBinContent(dataMaxBin)
          +h_realData.GetBinError(dataMaxBin);
        themax=std::max(themax,stack.GetMaximum());
        themax=std::max(themax,h_signal400.GetMaximum());
        themax=std::max(themax,h_signal250.GetMaximum());
        themax*=(1.0);
        stack.SetMaximum(themax);
        h_realData.SetMaximum(themax);
        h_signal400.SetMaximum(themax);
        h_signal250.SetMaximum(themax);

        std::ostringstream oss;
        if(signalScale!=1.0) oss<< "*" << signalScale;
        TCanvas canvas(obj_name.c_str(), obj_name.c_str());
        stack.Draw("hist");
        TLegend legend(0.7,0.5,0.95,0.85);
        legend.SetFillColor(4000);
        legend.AddEntry(&h_realData, "CMS Data, #sqrt{s}=8 TeV","lpe");
        legend.AddEntry(&h_TTJets_Hadronic,"ttbar (0l)","lf");
        legend.AddEntry(&h_TTJets_SemiLept,"ttbar (1l)","lf");
        legend.AddEntry(&h_TTJets_FullLept,"ttbar (2l)","lf");
        legend.AddEntry(&h_TTH,"ttH","lf");
        legend.AddEntry(&h_TTV,"ttV","lf");
        legend.AddEntry(&h_VH,"VH","lf");
        legend.AddEntry(&h_VV,"VV","lf");
        legend.AddEntry(&h_V,"V","lf");
        legend.AddEntry(&h_singleT, "single t", "lf");
        legend.AddEntry(&h_QCD,"QCD","lf");
        legend.AddEntry(&h_signal250,("Hbb(250)"+oss.str()).c_str(),"l");
        legend.AddEntry(&h_signal400,("Hbb(400)"+oss.str()).c_str(),"l");
        h_signal250.Draw("histsame");
        h_signal400.Draw("histsame");
        h_realData.Draw("e1psame");
        legend.Draw("same");
        canvas.Write();
        canvas.Print((argv[1]+std::string("_")+obj_name+".pdf").c_str());
        canvas.SetLogy(1);
        canvas.Print((argv[1]+std::string("_")+obj_name+"_log.pdf").c_str());
        canvas.SetLogy(0);
      }else if(QCD.at(0)->GetKey(obj_name.c_str(),1)->ReadObj()->IsA()==TH2D::Class()){
        TH2D h_QCD(GetTH2D(QCD,obj_name));
        TH2D h_TTJets_FullLept(GetTH2D(TTJets_FullLept,obj_name));
        TH2D h_TTJets_SemiLept(GetTH2D(TTJets_SemiLept,obj_name));
        TH2D h_TTJets_Hadronic(GetTH2D(TTJets_Hadronic,obj_name));
        TH2D h_TTH(GetTH2D(TTH,obj_name));
        TH2D h_TTV(GetTH2D(TTV,obj_name));
        TH2D h_VH(GetTH2D(VH,obj_name));
        TH2D h_VV(GetTH2D(VV,obj_name));
        TH2D h_V(GetTH2D(V,obj_name));
        TH2D h_singleT(GetTH2D(singleT,obj_name));
        TH2D h_realData(GetTH2D(realData,obj_name));
        TH2D h_signal400(GetTH2D(signal400,obj_name));
        TH2D h_signal250(GetTH2D(signal250,obj_name));

        TH2D h_backgroundMC=h_QCD+h_TTJets_FullLept+h_TTJets_SemiLept
          +h_TTJets_Hadronic+h_TTV+h_VH+h_VV+h_V+h_singleT;
        h_backgroundMC.SetFillColor(c_backgroundMC);
        h_realData.SetFillColor(c_realData);
        h_signal400.SetFillColor(c_signal400);
        h_signal250.SetFillColor(c_signal250);

        h_backgroundMC.SetLineColor(c_backgroundMC);
        h_realData.SetLineColor(c_realData);
        h_signal400.SetLineColor(c_signal400);
        h_signal250.SetLineColor(c_signal250);

        h_backgroundMC.SetMarkerColor(c_backgroundMC);
        h_realData.SetMarkerColor(c_realData);
        h_signal400.SetMarkerColor(c_signal400);
        h_signal250.SetMarkerColor(c_signal250);
        
        const double scaling=1.0;//h_realData.Integral("width")/(h_backgroundMC.Integral("width"));
        h_backgroundMC.Scale(scaling);
        h_signal400.Scale(scaling);
        h_signal250.Scale(scaling);

        h_signal400.Scale(signalScale);
        h_signal250.Scale(signalScale);

        h_backgroundMC.SetStats(false);
        h_realData.SetStats(false);
        h_signal400.SetStats(false);
        h_signal250.SetStats(false);

        std::ostringstream oss;
        if(signalScale!=1.0) oss << "*" << signalScale;

        TCanvas canvas(obj_name.c_str(), obj_name.c_str());
        h_signal400.Draw(GetScatArg(h_signal400).c_str());
        h_backgroundMC.Draw(GetScatArg(h_backgroundMC,"same").c_str());
        h_realData.Draw(GetScatArg(h_realData,"same").c_str());
        TLegend legend(0.8,0.85,1.0,1.0);
        legend.SetFillColor(4000);
        legend.AddEntry(&h_realData, "CMS Data, #sqrt{s}=8 TeV","f");
        legend.AddEntry(&h_backgroundMC,"SM","f");
        legend.AddEntry(&h_signal400,("Hbb(400)"+oss.str()).c_str(),"f");
        legend.Draw("same");
        canvas.Write();
        canvas.Print((argv[1]+std::string("_")+obj_name+".pdf").c_str());
      }
    }
    outFile.Close();
  }
}

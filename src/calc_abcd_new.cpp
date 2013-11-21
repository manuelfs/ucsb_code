#include <vector>
#include <iostream>
#include "stdint.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TMinuit.h"
#include "style.hpp"
#include "weights.hpp"
#include "math.hpp"
#include "minuit_functor.hpp"
#include "abcd_count.hpp"
#include "abcd_calculator.hpp"

const unsigned int max_draws(1000);
TRandom3 rng(3);

void GetFiles(std::vector<TFile*> &files, const std::vector<std::string> &names){
  for(std::vector<TFile*>::size_type file(0); file<files.size(); ++file){
    if(files.at(file)!=NULL && files.at(file)->IsOpen() && !files.at(file)->IsZombie()){
      files.at(file)->Close();
    }   
    delete files.at(file);
    files.at(file)=NULL;
  }
  files.clear();
  files.resize(names.size(), NULL);
  for(std::vector<std::string>::size_type name(0); name<names.size(); ++name){
    files.at(name)=new TFile(names.at(name).c_str(),"read");
  }
}

void GetTrees(std::vector<TTree*> &trees, const std::vector<TFile*> &files){
  for(std::vector<TTree*>::size_type tree(0); tree<trees.size(); ++tree){
    if(trees.at(tree)!=NULL){
      delete trees.at(tree);
      trees.at(tree)=NULL;
    }
  }
  trees.clear();
  trees.resize(files.size(), NULL);
  for(std::vector<TFile*>::size_type file(0); file<files.size(); ++file){
    if(files.at(file)!=NULL && files.at(file)->IsOpen() && !files.at(file)->IsZombie()){
      files.at(file)->GetObject("ABCD",trees.at(file));
    }
  }
}

void GetCounts(std::vector<ABCDCount>& counts, std::vector<TTree*>& trees,
	       const std::vector<double>& weights){
  counts.resize(trees.size());
  for(unsigned int i(0); i<trees.size(); ++i){
    if(trees.at(i)->GetEntries()>0){
      std::vector<double> a_count(4), b_count(4), c_count(4), d_count(4);
      trees.at(i)->SetBranchAddress("AWeighted_sbin1", &a_count.at(0));
      trees.at(i)->SetBranchAddress("AWeighted_sbin2", &a_count.at(1));
      trees.at(i)->SetBranchAddress("AWeighted_sbin3", &a_count.at(2));
      trees.at(i)->SetBranchAddress("AWeighted_sbin4", &a_count.at(3));
      trees.at(i)->SetBranchAddress("BWeighted_sbin1", &b_count.at(0));
      trees.at(i)->SetBranchAddress("BWeighted_sbin2", &b_count.at(1));
      trees.at(i)->SetBranchAddress("BWeighted_sbin3", &b_count.at(2));
      trees.at(i)->SetBranchAddress("BWeighted_sbin4", &b_count.at(3));
      trees.at(i)->SetBranchAddress("C2bWeighted_sbin1", &c_count.at(0));
      trees.at(i)->SetBranchAddress("C2bWeighted_sbin2", &c_count.at(1));
      trees.at(i)->SetBranchAddress("C2bWeighted_sbin3", &c_count.at(2));
      trees.at(i)->SetBranchAddress("C2bWeighted_sbin4", &c_count.at(3));
      trees.at(i)->SetBranchAddress("D2bWeighted_sbin1", &d_count.at(0));
      trees.at(i)->SetBranchAddress("D2bWeighted_sbin2", &d_count.at(1));
      trees.at(i)->SetBranchAddress("D2bWeighted_sbin3", &d_count.at(2));
      trees.at(i)->SetBranchAddress("D2bWeighted_sbin4", &d_count.at(3));
      trees.at(i)->GetEntry(0);
      counts.at(i).SetCounts(a_count, b_count, c_count, d_count);
      counts.at(i).SetWeight(weights.at(i));
    }
  }
}

void GetWeights(std::vector<double> &weights, const std::vector<std::string> &names,
                const WeightCalculator &weightCalc){
  weights.clear();
  weights.resize(names.size(), 0.0);
  for(std::vector<std::string>::size_type name(0); name<names.size(); ++name){
    weights.at(name)=weightCalc.GetWeight(names.at(name));
  }
}

void KillTrees(std::vector<TTree*> &trees){
  for(std::vector<TTree*>::size_type tree(0); tree<trees.size(); ++tree){
    if(trees.at(tree)!=NULL){
      delete trees.at(tree);
      trees.at(tree)=NULL;
    }
  }
}

void KillFiles(std::vector<TFile*> &files){
  for(std::vector<TFile*>::size_type file(0); file<files.size(); ++file){
    if(files.at(file)!=NULL && files.at(file)->IsOpen()){
      files.at(file)->Close();
    }
    delete files.at(file);
    files.at(file)=NULL;
  }
}

double GetRandom(const double kp1, const double weight){
  const double quantile(rng.Rndm());
  double left(0.0), right(DBL_MAX);
  double middle(left+0.5*(right-left));
  while(left<middle && middle<right){
    if(TMath::Gamma(kp1, middle)<=quantile){
      left=middle;
    }else{
      right=middle;
    }
    middle=left+0.5*(right-left);
  }
  return middle*weight;
}

void SetParameter(const unsigned int i, const unsigned int num_params, TMinuit& minuit){
  if(i==num_params-3){
    minuit.DefineParameter(i, "chi_b", 10.0, 1.0, 0.0, 0.0);
  }else if(i==num_params-2){
    minuit.DefineParameter(i, "chi_c", 10.0, 1.0, 0.0, 0.0);
  }else if(i==num_params-1){
    minuit.DefineParameter(i, "sig str", 0.5, 0.05, 0.0, 0.0);
  }else if(i<num_params-3){
    char name[128];
    switch(i%5){
    case 0:
      sprintf(name, "mu (%d)", i/5);
      break;
    case 1:
      sprintf(name, "As (%d)", i/5);
      break;
    case 2:
      sprintf(name, "Bs (%d)", i/5);
      break;
    case 3:
      sprintf(name, "Cs (%d)", i/5);
      break;
    case 4:
      sprintf(name, "Ds (%d)", i/5);
      break;
    }
    minuit.DefineParameter(i, name, 10.0, 1.0, 0.0, 0.0);
  }
}

int main(int argc, char *argv[]){
  SetStyle();
  bool use_mc(false);
  char opt(' ');
  while(( opt=getopt(argc, argv, "m") )!=-1){
    switch(opt){
    case 'm':
      use_mc=true;
      break;
    default:
      break;
    }
  }

  {TTree crap;}
  std::vector<std::string> observed_names(0), signal_names;
  signal_names.push_back("raw_plots_and_values/SMS-TChiHH_2b2b_2J_mChargino-400_mLSP-1_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1871_v71_SyncSkim.root");
  if(use_mc){
    observed_names.push_back("raw_plots_and_values/BJets_HT-1000ToInf_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1895_v71_SyncSkim.root");//6
    observed_names.push_back("raw_plots_and_values/BJets_HT-250To500_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1893_v71_SyncSkim.root");//7
    observed_names.push_back("raw_plots_and_values/BJets_HT-500To1000_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1894_v71_SyncSkim.root");//8
    observed_names.push_back("raw_plots_and_values/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71_SyncSkim.root");//9
    observed_names.push_back("raw_plots_and_values/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71_SyncSkim.root");//10
    observed_names.push_back("raw_plots_and_values/TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71_SyncSkim.root");//11
    observed_names.push_back("raw_plots_and_values/TTH_HToBB_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1855_v71_SyncSkim.root");//12
    observed_names.push_back("raw_plots_and_values/TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1857_v71_SyncSkim.root");//13
    observed_names.push_back("raw_plots_and_values/TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1856_v71_SyncSkim.root");//14
    observed_names.push_back("raw_plots_and_values/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1864_v71_SyncSkim.root");//15
    observed_names.push_back("raw_plots_and_values/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1865_v71_SyncSkim.root");//16
    observed_names.push_back("raw_plots_and_values/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1866_v71_SyncSkim.root");//17
    observed_names.push_back("raw_plots_and_values/T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1860_v71_SyncSkim.root");//18
    observed_names.push_back("raw_plots_and_values/T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1861_v71_SyncSkim.root");//19
    observed_names.push_back("raw_plots_and_values/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1862_v71_SyncSkim.root");//20
    observed_names.push_back("raw_plots_and_values/W2JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1877_v71_SyncSkim.root");//21
    observed_names.push_back("raw_plots_and_values/W3JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1878_v71_SyncSkim.root");//22
    observed_names.push_back("raw_plots_and_values/W4JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1879_v71_SyncSkim.root");//23
    observed_names.push_back("raw_plots_and_values/ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1887_v71_SyncSkim.root");//24
    observed_names.push_back("raw_plots_and_values/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1889_v71_SyncSkim.root");//25
    observed_names.push_back("raw_plots_and_values/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1888_v71_SyncSkim.root");//26
    observed_names.push_back("raw_plots_and_values/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1891_v71_SyncSkim.root");//27
    observed_names.push_back("raw_plots_and_values/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1890_v71_SyncSkim.root");//28
    observed_names.push_back("raw_plots_and_values/WH_WToLNu_HToBB_M-125_8TeV-powheg-herwigpp_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1858_v71_SyncSkim.root");//29
    observed_names.push_back("raw_plots_and_values/ZH_ZToBB_HToBB_M-125_8TeV-powheg-herwigpp_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1868_v71_SyncSkim.root");//30
    observed_names.push_back("raw_plots_and_values/WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1874_v71_SyncSkim.root");//31
    observed_names.push_back("raw_plots_and_values/ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1876_v71_SyncSkim.root");//32
    observed_names.push_back("raw_plots_and_values/WbbJetsToLNu_Massive_TuneZ2star_8TeV-madgraph-pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1859_v71_SyncSkim.root");//33
  }else{
    observed_names.push_back("raw_plots_and_values/MET_Run2012A-13Jul2012-v1_AOD_UCSB1852_v71_SyncSkim.root");//0
    observed_names.push_back("raw_plots_and_values/MET_Run2012B-13Jul2012-v1_AOD_UCSB1853_v71_SyncSkim.root");//1
    observed_names.push_back("raw_plots_and_values/MET_Run2012C-24Aug2012-v1_AOD_UCSB1854_v71_SyncSkim.root");//2
    observed_names.push_back("raw_plots_and_values/MET_Run2012C-PromptReco-v2_AOD_UCSB1867_v71_SyncSkim.root");//3
    observed_names.push_back("raw_plots_and_values/MET_Run2012D-PromptReco-v1_AOD_UCSB1869_v71_SyncSkim.root");//4
    observed_names.push_back("raw_plots_and_values/MET_Run2012D-PromptReco-v1_AOD_UCSB1870_v71_SyncSkim.root");//5
  }
  std::vector<TFile*> observed_files(0), signal_files(0);
  GetFiles(observed_files, observed_names);
  GetFiles(signal_files, signal_names);

  std::vector<TTree*> observed_trees(0), signal_trees(0);
  GetTrees(observed_trees, observed_files);
  GetTrees(signal_trees, signal_files);

  WeightCalculator weightCalc(19399.0);
  std::vector<double> observed_weights(0), signal_weights;
  GetWeights(observed_weights, observed_names, weightCalc);
  GetWeights(signal_weights, signal_names, weightCalc);

  std::vector<ABCDCount> observed_counts(0), signal_counts;
  GetCounts(observed_counts, observed_trees, observed_weights);
  GetCounts(signal_counts, signal_trees, signal_weights);

  KillTrees(observed_trees);
  KillTrees(signal_trees);
  KillFiles(observed_files);
  KillFiles(signal_files);

  ABCDCalculator abcd_calculator(observed_counts, signal_counts);
  const unsigned int num_params(abcd_calculator.GetNumberOfParameters());
  TMinuit minuit(num_params);
  minuit.SetPrintLevel(1);

  double strategy[1]={2.0};
  int useless(0);
  minuit.mnexcm("SET STR", strategy, 1, useless);
  MinuitFunctor<ABCDCalculator>::SetFunctor(&abcd_calculator);
  MinuitFunctor<ABCDCalculator>::SetNumParams(num_params);
  minuit.SetFCN(MinuitFunctor<ABCDCalculator>::Function);
  minuit.SetMaxIterations(std::numeric_limits<int32_t>::max());
  std::cout << "before" << std::endl;
  for(unsigned int i(0); i<num_params; ++i){
    SetParameter(i, num_params, minuit);
  }
  std::cout << "middle" << std::endl;
  minuit.Migrad();
  minuit.Migrad();
  minuit.Migrad();
  minuit.mnmnos();
  minuit.mnmnos();
  minuit.mnmnos();
  std::cout << "after" << std::endl;
}

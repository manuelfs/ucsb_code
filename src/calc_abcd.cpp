#include <cmath>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <string>
#include <unistd.h>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "style.hpp"
#include "weights.hpp"

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

void GetValues(std::vector<double> &A, std::vector<double> &B, std::vector<double> &C3,
               std::vector<double> &D3, std::vector<double> &C2, std::vector<double> &D2,
               const std::vector<TTree*> &trees, const std::string &str=""){
  A.clear();
  B.clear();
  C3.clear();
  D3.clear();
  C2.clear();
  D2.clear();
  A.resize(trees.size(),0.0);
  B.resize(trees.size(),0.0);
  C3.resize(trees.size(),0.0);
  D3.resize(trees.size(),0.0);
  C2.resize(trees.size(),0.0);
  D2.resize(trees.size(),0.0);
  for(std::vector<TTree*>::size_type tree(0); tree<trees.size(); ++tree){
    if(trees.at(tree)!=NULL && trees.at(tree)->GetEntries()>0){
      trees.at(tree)->SetBranchAddress(("AWeighted"+str).c_str(),&A.at(tree));
      trees.at(tree)->SetBranchAddress(("BWeighted"+str).c_str(),&B.at(tree));
      trees.at(tree)->SetBranchAddress(("C3bWeighted"+str).c_str(),&C3.at(tree));
      trees.at(tree)->SetBranchAddress(("D3bWeighted"+str).c_str(),&D3.at(tree));
      trees.at(tree)->SetBranchAddress(("C2bWeighted"+str).c_str(),&C2.at(tree));
      trees.at(tree)->SetBranchAddress(("D2bWeighted"+str).c_str(),&D2.at(tree));
      
      trees.at(tree)->GetEntry(0);
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

void GetUncertaintiesOld(double &up, double &down, const double center,
                         const std::vector<double> vals){
  if(vals.size()!=0){
    const double frac(TMath::Erf(1.0/sqrt(2.0)));
    const std::vector<double>::size_type delta(static_cast<int>(ceil(frac*(vals.size()-1))));
    double min_diff(DBL_MAX);
    double left(0.0), right(0.0);
    for(std::vector<double>::size_type val(0); val+delta<vals.size(); ++val){
      std::vector<double>::size_type val2(val+delta);
      const double diff(vals.at(val2)-vals.at(val));
      if(diff<min_diff){
        min_diff=diff;
        right=vals.at(val2);
        left=vals.at(val);
      }
    }
    if(center<left || center>right){
      std::cout << left << " " << center << " " << right << std::endl;
    }
    up=right-center;
    down=center-left;
    if(up<0.0) up=0.0;
    if(down<0.0) down=0.0;
  }else{
    up=0.0;
    down=0.0;
  }
}

void GetUncertainties(double &up, double &down, const double center,
                      const std::vector<double> vals){
  unsigned int best_pos(-1);
  double best_diff(DBL_MAX);
  for(unsigned int pos(0); pos<vals.size(); ++pos){
    const double diff(fabs(center-vals.at(pos)));
    if(diff<best_diff){
      best_diff=diff;
      best_pos=pos;
    }
  }
  const double frac(TMath::Erf(1.0/sqrt(2.0)));
  const std::vector<double>::size_type delta(static_cast<int>(0.5*ceil(frac*(vals.size()-1))));
  unsigned int upper(best_pos+delta), lower(best_pos-delta);
  if(lower>best_pos) lower=best_pos;
  if(upper>=vals.size()) upper=vals.size()-1;
  up=vals.at(upper)-vals.at(best_pos);
  down=vals.at(best_pos)-vals.at(lower);
}

void GetSummedVals(double &up, double &down, double &center,
                   const std::vector<double> &val,
                   const std::vector<double> &weights,
                   const std::vector<double>::size_type &low,
                   const std::vector<double>::size_type &high,
                   const bool verbose=false){
  center=0.0;
  for(std::vector<double>::size_type sample(low); sample<high; ++sample){
    center+=val.at(sample);
  }
  std::vector<double> draws(0);
  for(unsigned int draw(0); draw<max_draws; ++draw){
    double thisVal(0.0);
    for(std::vector<double>::size_type sample(low); sample<high; ++sample){
      const double blah(GetRandom(val.at(sample)/weights.at(sample), weights.at(sample)));
      if(verbose){
        std::cout << sample << " " << val.at(sample) << " " <<  weights.at(sample) << " " << blah << std::endl;
      }
      thisVal+=blah;
    }
    if(verbose) std::cout << thisVal << std::endl;
    draws.push_back(thisVal);
  }
  std::sort(draws.begin(), draws.end());
  GetUncertainties(up, down, center, draws);
}

void GetSummedKappas(double &up, double &down, double &center, const std::vector<double> &A,
                     const std::vector<double> &B, const std::vector<double> &C,
                     const std::vector<double> &D, const std::vector<double> &weights,
                     const std::vector<double>::size_type &low,
                     const std::vector<double>::size_type &high){
  std::vector<double> draws(0);
  for(unsigned int draw(0); draw<max_draws; ++draw){
    double sumA(0.0), sumB(0.0), sumC(0.0), sumD(0.0);
    for(std::vector<double>::size_type sample(low); sample<A.size() && sample<high; ++sample){
      sumA+=GetRandom(A.at(sample)/weights.at(sample), weights.at(sample));
      sumB+=GetRandom(B.at(sample)/weights.at(sample), weights.at(sample));
      sumC+=GetRandom(C.at(sample)/weights.at(sample), weights.at(sample));
      sumD+=GetRandom(D.at(sample)/weights.at(sample), weights.at(sample));
    }
    if(sumB*sumC>0.0){
      draws.push_back(sumA*sumD/(sumB*sumC));
    }else{
      draws.push_back(DBL_MAX);
    }
  }

  double sumA(0.0), sumB(0.0), sumC(0.0), sumD(0.0);
  for(std::vector<double>::size_type sample(low); sample<high; ++sample){
    sumA+=A.at(sample);
    sumB+=B.at(sample);
    sumC+=C.at(sample);
    sumD+=D.at(sample);
  }
  if(sumB*sumC>0.0){
    center=sumA*sumD/(sumB*sumC);
  }else{
    center=DBL_MAX;
  }
  std::sort(draws.begin(), draws.end());
  GetUncertainties(up, down, center, draws);
}

void GetSummedAPred(double &up, double &down, double &center,
                    const std::vector<double> &B, const std::vector<double> &C,
                    const std::vector<double> &D, const std::vector<double> &weights,
                    const std::vector<double>::size_type &low,
                    const std::vector<double>::size_type &high,
                    const bool verbose=false){
  std::vector<double> draws(0);
  for(unsigned int draw(0); draw<max_draws; ++draw){
    double sumB(0.0), sumC(0.0), sumD(0.0);
    for(std::vector<double>::size_type sample(low); sample<B.size() && sample<high; ++sample){
      sumB+=GetRandom(B.at(sample)/weights.at(sample), weights.at(sample));
      sumC+=GetRandom(C.at(sample)/weights.at(sample), weights.at(sample));
      sumD+=GetRandom(D.at(sample)/weights.at(sample), weights.at(sample));
    }
    if(verbose) std::cout << sumB << " " << sumC << " " << sumD << std::endl;
    if(sumD>0.0){
      draws.push_back(sumB*sumC/sumD);
    }else{
      draws.push_back(DBL_MAX);
    }
  }

  double sumB(0.0), sumC(0.0), sumD(0.0);
  for(std::vector<double>::size_type sample(low); sample<high; ++sample){
    sumB+=B.at(sample);
    sumC+=C.at(sample);
    sumD+=D.at(sample);
  }
  if(sumD>0.0){
    center=sumB*sumC/sumD;
  }else{
    center=DBL_MAX;
  }
  std::sort(draws.begin(), draws.end());
  GetUncertainties(up, down, center, draws);
}

void PrintLine(std::string name, const double A, const double Aup, const double Adown, const double B, const double Bup, const double Bdown, const double C3, const double C3up, const double C3down, const double D3, const double D3up, const double D3down, const double C2, const double C2up, const double C2down, const double D2, const double D2up, const double D2down, const double kappa3, const double kappa3up, const double kappa3down, const double kappa2, const double kappa2up, const double kappa2down, const double pred23, const double predup23, const double preddown23, const double pred24, const double predup24, const double preddown24, const double pred34, const double predup34, const double preddown34){
  std::cout << name << " & $"
            << std::setprecision(3) << A << "_{-"
            << std::setprecision(2) << Adown << "}^{+" << Aup << "}$ & $"
            << std::setprecision(3) << B << "_{-"
            << std::setprecision(2) << Bdown << "}^{+" << Bup << "}$ & $"
            << std::setprecision(3) << C3 << "_{-"
            << std::setprecision(2) << C3down << "}^{+" << C3up << "}$ & $"
            << std::setprecision(3) << D3 << "_{-"
            << std::setprecision(2) << D3down << "}^{+" << D3up << "}$ & $"
            << std::setprecision(3) << C2 << "_{-"
            << std::setprecision(2) << C2down << "}^{+" << C2up << "}$ & $"
            << std::setprecision(3) << D2 << "_{-"
            << std::setprecision(2) << D2down << "}^{+" << D2up << "}$ & $"
            << std::setprecision(3) << kappa3 << "_{-"
            << std::setprecision(2) << kappa3down << "}^{+" << kappa3up << "}$ & $"
            << std::setprecision(3) << kappa2 << "_{-"
            << std::setprecision(2) << kappa2down << "}^{+" << kappa2up << "}$ & $"
            << std::setprecision(3) << pred23 << "_{-"
            << std::setprecision(2) << preddown23 << "}^{+" << predup23 << "}$ & $"
            << std::setprecision(3) << pred24 << "_{-"
            << std::setprecision(2) << preddown24 << "}^{+" << predup24 << "}$ & $"
            << std::setprecision(3) << pred34 << "_{-"
            << std::setprecision(2) << preddown34 << "}^{+" << predup34 << "}$"
            << std::endl;
}

int main(int argc, char *argv[]){
  SetStyle();
  bool plot_only(false), mc_plot(false);
  char opt(' ');
  while(( opt=getopt(argc, argv, "pm") )!=-1){
    switch(opt){
    case 'p':
      plot_only=true;
      break;
    case 'm':
      mc_plot=true;
      break;
    default:
      std::cerr << "Error in " << argv[0] << ": '" << opt
                << "' is not a valid option." << std::endl;
    }
  }

  {TTree crap;}
  std::vector<std::string> names(0);
  names.push_back("raw_plots_and_values/MET_Run2012A-13Jul2012-v1_AOD_UCSB1852_v71_SyncSkim.root");//0
  names.push_back("raw_plots_and_values/MET_Run2012B-13Jul2012-v1_AOD_UCSB1853_v71_SyncSkim.root");//1
  names.push_back("raw_plots_and_values/MET_Run2012C-24Aug2012-v1_AOD_UCSB1854_v71_SyncSkim.root");//2
  names.push_back("raw_plots_and_values/MET_Run2012C-PromptReco-v2_AOD_UCSB1867_v71_SyncSkim.root");//3
  names.push_back("raw_plots_and_values/MET_Run2012D-PromptReco-v1_AOD_UCSB1869_v71_SyncSkim.root");//4
  names.push_back("raw_plots_and_values/MET_Run2012D-PromptReco-v1_AOD_UCSB1870_v71_SyncSkim.root");//5
  names.push_back("raw_plots_and_values/BJets_HT-1000ToInf_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1895_v71_SyncSkim.root");//6
  names.push_back("raw_plots_and_values/BJets_HT-250To500_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1893_v71_SyncSkim.root");//7
  names.push_back("raw_plots_and_values/BJets_HT-500To1000_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1894_v71_SyncSkim.root");//8
  names.push_back("raw_plots_and_values/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71_SyncSkim.root");//9
  names.push_back("raw_plots_and_values/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71_SyncSkim.root");//10
  names.push_back("raw_plots_and_values/TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71_SyncSkim.root");//11
  names.push_back("raw_plots_and_values/TTH_HToBB_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1855_v71_SyncSkim.root");//12
  names.push_back("raw_plots_and_values/TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1857_v71_SyncSkim.root");//13
  names.push_back("raw_plots_and_values/TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1856_v71_SyncSkim.root");//14
  names.push_back("raw_plots_and_values/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1864_v71_SyncSkim.root");//15
  names.push_back("raw_plots_and_values/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1865_v71_SyncSkim.root");//16
  names.push_back("raw_plots_and_values/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1866_v71_SyncSkim.root");//17
  names.push_back("raw_plots_and_values/T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1860_v71_SyncSkim.root");//18
  names.push_back("raw_plots_and_values/T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1861_v71_SyncSkim.root");//19
  names.push_back("raw_plots_and_values/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1862_v71_SyncSkim.root");//20
  names.push_back("raw_plots_and_values/W2JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1877_v71_SyncSkim.root");//21
  names.push_back("raw_plots_and_values/W3JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1878_v71_SyncSkim.root");//22
  names.push_back("raw_plots_and_values/W4JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1879_v71_SyncSkim.root");//23
  names.push_back("raw_plots_and_values/ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1887_v71_SyncSkim.root");//24
  names.push_back("raw_plots_and_values/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1889_v71_SyncSkim.root");//25
  names.push_back("raw_plots_and_values/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1888_v71_SyncSkim.root");//26
  names.push_back("raw_plots_and_values/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1891_v71_SyncSkim.root");//27
  names.push_back("raw_plots_and_values/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1890_v71_SyncSkim.root");//28
  names.push_back("raw_plots_and_values/WH_WToLNu_HToBB_M-125_8TeV-powheg-herwigpp_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1858_v71_SyncSkim.root");//29
  names.push_back("raw_plots_and_values/ZH_ZToBB_HToBB_M-125_8TeV-powheg-herwigpp_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1868_v71_SyncSkim.root");//30
  names.push_back("raw_plots_and_values/WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1874_v71_SyncSkim.root");//31
  names.push_back("raw_plots_and_values/ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1876_v71_SyncSkim.root");//32
  names.push_back("raw_plots_and_values/WbbJetsToLNu_Massive_TuneZ2star_8TeV-madgraph-pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1859_v71_SyncSkim.root");//33

  std::vector<TFile*> files(0);
  GetFiles(files, names);

  std::vector<TTree*> trees(0);
  GetTrees(trees, files);

  std::vector<double> A(0), B(0), C3(0), D3(0), C2(0), D2(0);
  GetValues(A, B, C3, D3, C2, D2, trees, "");

  std::vector<double> A_sbin1(0), B_sbin1(0), C3_sbin1(0), D3_sbin1(0), C2_sbin1(0), D2_sbin1(0);
  GetValues(A_sbin1, B_sbin1, C3_sbin1, D3_sbin1, C2_sbin1, D2_sbin1, trees, "_sbin1");
  std::vector<double> A_sbin2(0), B_sbin2(0), C3_sbin2(0), D3_sbin2(0), C2_sbin2(0), D2_sbin2(0);
  GetValues(A_sbin2, B_sbin2, C3_sbin2, D3_sbin2, C2_sbin2, D2_sbin2, trees, "_sbin2");
  std::vector<double> A_sbin3(0), B_sbin3(0), C3_sbin3(0), D3_sbin3(0), C2_sbin3(0), D2_sbin3(0);
  GetValues(A_sbin3, B_sbin3, C3_sbin3, D3_sbin3, C2_sbin3, D2_sbin3, trees, "_sbin3");
  std::vector<double> A_sbin4(0), B_sbin4(0), C3_sbin4(0), D3_sbin4(0), C2_sbin4(0), D2_sbin4(0);
  GetValues(A_sbin4, B_sbin4, C3_sbin4, D3_sbin4, C2_sbin4, D2_sbin4, trees, "_sbin4");

  WeightCalculator weightCalc(19399.0);
  std::vector<double> weights(0);
  GetWeights(weights, names, weightCalc);

  std::vector<std::string> tex(0);
  std::vector<unsigned int> upper(0), lower(0);
  if(!(plot_only && mc_plot)){
    tex.push_back("Data"); upper.push_back(6); lower.push_back(0);
  }
  if(!plot_only){
    tex.push_back("QCD"); upper.push_back(9); lower.push_back(6);
    tex.push_back("$t\\overline{t}$ (2l)"); upper.push_back(10); lower.push_back(9);
    tex.push_back("$t\\overline{t}$ (1l)"); upper.push_back(11); lower.push_back(10);
    tex.push_back("$t\\overline{t}$ (0l)"); upper.push_back(12); lower.push_back(11);
    tex.push_back("t\\overline{t}H$"); upper.push_back(13); lower.push_back(12);
    tex.push_back("$t\\overline{t}V$"); upper.push_back(15); lower.push_back(13);
    tex.push_back("$t$"); upper.push_back(21); lower.push_back(15);
    tex.push_back("$V$"); upper.push_back(29); lower.push_back(21);
    tex.push_back("$VH$"); upper.push_back(31); lower.push_back(29);
    tex.push_back("$VV$"); upper.push_back(33); lower.push_back(31);
  }
  if(mc_plot || !plot_only){
    tex.push_back("SM total"); upper.push_back(33); lower.push_back(9);
  }

  std::vector<double> Aval(tex.size()), Aup(tex.size()), Adown(tex.size()),
    Bval(tex.size()), Bup(tex.size()), Bdown(tex.size()),
    C3val(tex.size()), C3up(tex.size()), C3down(tex.size()),
    D3val(tex.size()), D3up(tex.size()), D3down(tex.size()),
    C2val(tex.size()), C2up(tex.size()), C2down(tex.size()),
    D2val(tex.size()), D2up(tex.size()), D2down(tex.size()),
    kappa3val(tex.size()), kappa3up(tex.size()), kappa3down(tex.size()),
    kappa2val(tex.size()), kappa2up(tex.size()), kappa2down(tex.size());

  std::vector<double> Aval_sbin1(tex.size()), Aup_sbin1(tex.size()), Adown_sbin1(tex.size()),
    Bval_sbin1(tex.size()), Bup_sbin1(tex.size()), Bdown_sbin1(tex.size()),
    C3val_sbin1(tex.size()), C3up_sbin1(tex.size()), C3down_sbin1(tex.size()),
    D3val_sbin1(tex.size()), D3up_sbin1(tex.size()), D3down_sbin1(tex.size()),
    C2val_sbin1(tex.size()), C2up_sbin1(tex.size()), C2down_sbin1(tex.size()),
    D2val_sbin1(tex.size()), D2up_sbin1(tex.size()), D2down_sbin1(tex.size());

  std::vector<double> Aval_sbin2(tex.size()), Aup_sbin2(tex.size()), Adown_sbin2(tex.size()),
    Bval_sbin2(tex.size()), Bup_sbin2(tex.size()), Bdown_sbin2(tex.size()),
    C3val_sbin2(tex.size()), C3up_sbin2(tex.size()), C3down_sbin2(tex.size()),
    D3val_sbin2(tex.size()), D3up_sbin2(tex.size()), D3down_sbin2(tex.size()),
    C2val_sbin2(tex.size()), C2up_sbin2(tex.size()), C2down_sbin2(tex.size()),
    D2val_sbin2(tex.size()), D2up_sbin2(tex.size()), D2down_sbin2(tex.size());

  std::vector<double> Aval_sbin3(tex.size()), Aup_sbin3(tex.size()), Adown_sbin3(tex.size()),
    Bval_sbin3(tex.size()), Bup_sbin3(tex.size()), Bdown_sbin3(tex.size()),
    C3val_sbin3(tex.size()), C3up_sbin3(tex.size()), C3down_sbin3(tex.size()),
    D3val_sbin3(tex.size()), D3up_sbin3(tex.size()), D3down_sbin3(tex.size()),
    C2val_sbin3(tex.size()), C2up_sbin3(tex.size()), C2down_sbin3(tex.size()),
    D2val_sbin3(tex.size()), D2up_sbin3(tex.size()), D2down_sbin3(tex.size());

  std::vector<double> Aval_sbin4(tex.size()), Aup_sbin4(tex.size()), Adown_sbin4(tex.size()),
    Bval_sbin4(tex.size()), Bup_sbin4(tex.size()), Bdown_sbin4(tex.size()),
    C3val_sbin4(tex.size()), C3up_sbin4(tex.size()), C3down_sbin4(tex.size()),
    D3val_sbin4(tex.size()), D3up_sbin4(tex.size()), D3down_sbin4(tex.size()),
    C2val_sbin4(tex.size()), C2up_sbin4(tex.size()), C2down_sbin4(tex.size()),
    D2val_sbin4(tex.size()), D2up_sbin4(tex.size()), D2down_sbin4(tex.size());

  std::vector<double> pred23(tex.size()), predup23(tex.size()), preddown23(tex.size());
  std::vector<double> pred24(tex.size()), predup24(tex.size()), preddown24(tex.size());
  std::vector<double> pred34(tex.size()), predup34(tex.size()), preddown34(tex.size());
  std::vector<double> pred32(tex.size()), predup32(tex.size()), preddown32(tex.size());
  std::vector<double> pred42(tex.size()), predup42(tex.size()), preddown42(tex.size());
  std::vector<double> pred43(tex.size()), predup43(tex.size()), preddown43(tex.size());

  std::vector<double> pred23_sbin1(tex.size()), predup23_sbin1(tex.size()), preddown23_sbin1(tex.size());
  std::vector<double> pred24_sbin1(tex.size()), predup24_sbin1(tex.size()), preddown24_sbin1(tex.size());
  std::vector<double> pred34_sbin1(tex.size()), predup34_sbin1(tex.size()), preddown34_sbin1(tex.size());
  std::vector<double> pred32_sbin1(tex.size()), predup32_sbin1(tex.size()), preddown32_sbin1(tex.size());
  std::vector<double> pred42_sbin1(tex.size()), predup42_sbin1(tex.size()), preddown42_sbin1(tex.size());
  std::vector<double> pred43_sbin1(tex.size()), predup43_sbin1(tex.size()), preddown43_sbin1(tex.size());

  std::vector<double> pred23_sbin2(tex.size()), predup23_sbin2(tex.size()), preddown23_sbin2(tex.size());
  std::vector<double> pred24_sbin2(tex.size()), predup24_sbin2(tex.size()), preddown24_sbin2(tex.size());
  std::vector<double> pred34_sbin2(tex.size()), predup34_sbin2(tex.size()), preddown34_sbin2(tex.size());
  std::vector<double> pred32_sbin2(tex.size()), predup32_sbin2(tex.size()), preddown32_sbin2(tex.size());
  std::vector<double> pred42_sbin2(tex.size()), predup42_sbin2(tex.size()), preddown42_sbin2(tex.size());
  std::vector<double> pred43_sbin2(tex.size()), predup43_sbin2(tex.size()), preddown43_sbin2(tex.size());

  std::vector<double> pred23_sbin3(tex.size()), predup23_sbin3(tex.size()), preddown23_sbin3(tex.size());
  std::vector<double> pred24_sbin3(tex.size()), predup24_sbin3(tex.size()), preddown24_sbin3(tex.size());
  std::vector<double> pred34_sbin3(tex.size()), predup34_sbin3(tex.size()), preddown34_sbin3(tex.size());
  std::vector<double> pred32_sbin3(tex.size()), predup32_sbin3(tex.size()), preddown32_sbin3(tex.size());
  std::vector<double> pred42_sbin3(tex.size()), predup42_sbin3(tex.size()), preddown42_sbin3(tex.size());
  std::vector<double> pred43_sbin3(tex.size()), predup43_sbin3(tex.size()), preddown43_sbin3(tex.size());

  std::vector<double> pred23_sbin4(tex.size()), predup23_sbin4(tex.size()), preddown23_sbin4(tex.size());
  std::vector<double> pred24_sbin4(tex.size()), predup24_sbin4(tex.size()), preddown24_sbin4(tex.size());
  std::vector<double> pred34_sbin4(tex.size()), predup34_sbin4(tex.size()), preddown34_sbin4(tex.size());
  std::vector<double> pred32_sbin4(tex.size()), predup32_sbin4(tex.size()), preddown32_sbin4(tex.size());
  std::vector<double> pred42_sbin4(tex.size()), predup42_sbin4(tex.size()), preddown42_sbin4(tex.size());
  std::vector<double> pred43_sbin4(tex.size()), predup43_sbin4(tex.size()), preddown43_sbin4(tex.size());

  for(unsigned int i(0); i<Aval.size(); ++i){
    GetSummedVals(Aup.at(i), Adown.at(i), Aval.at(i), A,
                  weights, lower.at(i), upper.at(i));
    GetSummedVals(Bup.at(i), Bdown.at(i), Bval.at(i), B,
                  weights, lower.at(i), upper.at(i));
    GetSummedVals(C3up.at(i), C3down.at(i), C3val.at(i), C3,
                  weights, lower.at(i), upper.at(i));
    GetSummedVals(D3up.at(i), D3down.at(i), D3val.at(i), D3,
                  weights, lower.at(i), upper.at(i));
    GetSummedVals(C2up.at(i), C2down.at(i), C2val.at(i), C2,
                  weights, lower.at(i), upper.at(i));
    GetSummedVals(D2up.at(i), D2down.at(i), D2val.at(i), D2,
                  weights, lower.at(i), upper.at(i));

    if(plot_only){
      GetSummedVals(Aup_sbin1.at(i), Adown_sbin1.at(i), Aval_sbin1.at(i), A_sbin1,
                    weights, lower.at(i), upper.at(i));
      GetSummedVals(Bup_sbin1.at(i), Bdown_sbin1.at(i), Bval_sbin1.at(i), B_sbin1,
                    weights, lower.at(i), upper.at(i));
      GetSummedVals(C3up_sbin1.at(i), C3down_sbin1.at(i), C3val_sbin1.at(i), C3_sbin1,
                    weights, lower.at(i), upper.at(i));
      GetSummedVals(D3up_sbin1.at(i), D3down_sbin1.at(i), D3val_sbin1.at(i), D3_sbin1,
                    weights, lower.at(i), upper.at(i));
      GetSummedVals(C2up_sbin1.at(i), C2down_sbin1.at(i), C2val_sbin1.at(i), C2_sbin1,
                    weights, lower.at(i), upper.at(i));
      GetSummedVals(D2up_sbin1.at(i), D2down_sbin1.at(i), D2val_sbin1.at(i), D2_sbin1,
                    weights, lower.at(i), upper.at(i));
      
      GetSummedVals(Aup_sbin2.at(i), Adown_sbin2.at(i), Aval_sbin2.at(i), A_sbin2,
                    weights, lower.at(i), upper.at(i));
      GetSummedVals(Bup_sbin2.at(i), Bdown_sbin2.at(i), Bval_sbin2.at(i), B_sbin2,
                    weights, lower.at(i), upper.at(i));
      GetSummedVals(C3up_sbin2.at(i), C3down_sbin2.at(i), C3val_sbin2.at(i), C3_sbin2,
                    weights, lower.at(i), upper.at(i));
      GetSummedVals(D3up_sbin2.at(i), D3down_sbin2.at(i), D3val_sbin2.at(i), D3_sbin2,
                    weights, lower.at(i), upper.at(i));
      GetSummedVals(C2up_sbin2.at(i), C2down_sbin2.at(i), C2val_sbin2.at(i), C2_sbin2,
                    weights, lower.at(i), upper.at(i));
      GetSummedVals(D2up_sbin2.at(i), D2down_sbin2.at(i), D2val_sbin2.at(i), D2_sbin2,
                    weights, lower.at(i), upper.at(i));
      
      GetSummedVals(Aup_sbin3.at(i), Adown_sbin3.at(i), Aval_sbin3.at(i), A_sbin3,
                    weights, lower.at(i), upper.at(i));
      GetSummedVals(Bup_sbin3.at(i), Bdown_sbin3.at(i), Bval_sbin3.at(i), B_sbin3,
                    weights, lower.at(i), upper.at(i));
      GetSummedVals(C3up_sbin3.at(i), C3down_sbin3.at(i), C3val_sbin3.at(i), C3_sbin3,
                    weights, lower.at(i), upper.at(i));
      GetSummedVals(D3up_sbin3.at(i), D3down_sbin3.at(i), D3val_sbin3.at(i), D3_sbin3,
                    weights, lower.at(i), upper.at(i));
      GetSummedVals(C2up_sbin3.at(i), C2down_sbin3.at(i), C2val_sbin3.at(i), C2_sbin3,
                    weights, lower.at(i), upper.at(i));
      GetSummedVals(D2up_sbin3.at(i), D2down_sbin3.at(i), D2val_sbin3.at(i), D2_sbin3,
                    weights, lower.at(i), upper.at(i));
      
      GetSummedVals(Aup_sbin4.at(i), Adown_sbin4.at(i), Aval_sbin4.at(i), A_sbin4,
                    weights, lower.at(i), upper.at(i));
      GetSummedVals(Bup_sbin4.at(i), Bdown_sbin4.at(i), Bval_sbin4.at(i), B_sbin4,
                    weights, lower.at(i), upper.at(i));
      GetSummedVals(C3up_sbin4.at(i), C3down_sbin4.at(i), C3val_sbin4.at(i), C3_sbin4,
                    weights, lower.at(i), upper.at(i));
      GetSummedVals(D3up_sbin4.at(i), D3down_sbin4.at(i), D3val_sbin4.at(i), D3_sbin4,
                    weights, lower.at(i), upper.at(i));
      GetSummedVals(C2up_sbin4.at(i), C2down_sbin4.at(i), C2val_sbin4.at(i), C2_sbin4,
                    weights, lower.at(i), upper.at(i));
      GetSummedVals(D2up_sbin4.at(i), D2down_sbin4.at(i), D2val_sbin4.at(i), D2_sbin4,
                    weights, lower.at(i), upper.at(i));
    }

    GetSummedKappas(kappa3up.at(i), kappa3down.at(i), kappa3val.at(i), A, B, C3, D3,
                    weights, lower.at(i), upper.at(i));
    GetSummedKappas(kappa2up.at(i), kappa2down.at(i), kappa2val.at(i), A, B, C2, D2,
                    weights, lower.at(i), upper.at(i));

    GetSummedAPred(predup23.at(i), preddown23.at(i), pred23.at(i), D3, C2, D2, weights,
                   lower.at(i), upper.at(i));
    GetSummedAPred(predup24.at(i), preddown24.at(i), pred24.at(i), B, C2, D2, weights,
                   lower.at(i), upper.at(i));
    GetSummedAPred(predup34.at(i), preddown34.at(i), pred34.at(i), B, C3, D3, weights,
                   lower.at(i), upper.at(i));
    if(plot_only){
      GetSummedAPred(predup32.at(i), preddown32.at(i), pred32.at(i), D2, C3, D3, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup42.at(i), preddown42.at(i), pred42.at(i), D2, A, B, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup43.at(i), preddown43.at(i), pred43.at(i), D3, A, B, weights,
                     lower.at(i), upper.at(i));
    }

    if(plot_only){
      GetSummedAPred(predup23_sbin1.at(i), preddown23_sbin1.at(i), pred23_sbin1.at(i), D3_sbin1, C2_sbin1, D2_sbin1, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup24_sbin1.at(i), preddown24_sbin1.at(i), pred24_sbin1.at(i), B_sbin1, C2_sbin1, D2_sbin1, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup34_sbin1.at(i), preddown34_sbin1.at(i), pred34_sbin1.at(i), B_sbin1, C3_sbin1, D3_sbin1, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup32_sbin1.at(i), preddown32_sbin1.at(i), pred32_sbin1.at(i), D2_sbin1, C3_sbin1, D3_sbin1, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup42_sbin1.at(i), preddown42_sbin1.at(i), pred42_sbin1.at(i), D2_sbin1, A_sbin1, B_sbin1, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup43_sbin1.at(i), preddown43_sbin1.at(i), pred43_sbin1.at(i), D3_sbin1, A_sbin1, B_sbin1, weights,
                     lower.at(i), upper.at(i));
      
      GetSummedAPred(predup23_sbin2.at(i), preddown23_sbin2.at(i), pred23_sbin2.at(i), D3_sbin2, C2_sbin2, D2_sbin2, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup24_sbin2.at(i), preddown24_sbin2.at(i), pred24_sbin2.at(i), B_sbin2, C2_sbin2, D2_sbin2, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup34_sbin2.at(i), preddown34_sbin2.at(i), pred34_sbin2.at(i), B_sbin2, C3_sbin2, D3_sbin2, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup32_sbin2.at(i), preddown32_sbin2.at(i), pred32_sbin2.at(i), D2_sbin2, C3_sbin2, D3_sbin2, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup42_sbin2.at(i), preddown42_sbin2.at(i), pred42_sbin2.at(i), D2_sbin2, A_sbin2, B_sbin2, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup43_sbin2.at(i), preddown43_sbin2.at(i), pred43_sbin2.at(i), D3_sbin2, A_sbin2, B_sbin2, weights,
                     lower.at(i), upper.at(i));
      
      GetSummedAPred(predup23_sbin3.at(i), preddown23_sbin3.at(i), pred23_sbin3.at(i), D3_sbin3, C2_sbin3, D2_sbin3, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup24_sbin3.at(i), preddown24_sbin3.at(i), pred24_sbin3.at(i), B_sbin3, C2_sbin3, D2_sbin3, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup34_sbin3.at(i), preddown34_sbin3.at(i), pred34_sbin3.at(i), B_sbin3, C3_sbin3, D3_sbin3, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup32_sbin3.at(i), preddown32_sbin3.at(i), pred32_sbin3.at(i), D2_sbin3, C3_sbin3, D3_sbin3, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup42_sbin3.at(i), preddown42_sbin3.at(i), pred42_sbin3.at(i), D2_sbin3, A_sbin3, B_sbin3, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup43_sbin3.at(i), preddown43_sbin3.at(i), pred43_sbin3.at(i), D3_sbin3, A_sbin3, B_sbin3, weights,
                     lower.at(i), upper.at(i));
      
      GetSummedAPred(predup23_sbin4.at(i), preddown23_sbin4.at(i), pred23_sbin4.at(i), D3_sbin4, C2_sbin4, D2_sbin4, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup24_sbin4.at(i), preddown24_sbin4.at(i), pred24_sbin4.at(i), B_sbin4, C2_sbin4, D2_sbin4, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup34_sbin4.at(i), preddown34_sbin4.at(i), pred34_sbin4.at(i), B_sbin4, C3_sbin4, D3_sbin4, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup32_sbin4.at(i), preddown32_sbin4.at(i), pred32_sbin4.at(i), D2_sbin4, C3_sbin4, D3_sbin4, weights,
                     lower.at(i), upper.at(i), true);
      GetSummedAPred(predup42_sbin4.at(i), preddown42_sbin4.at(i), pred42_sbin4.at(i), D2_sbin4, A_sbin4, B_sbin4, weights,
                     lower.at(i), upper.at(i));
      GetSummedAPred(predup43_sbin4.at(i), preddown43_sbin4.at(i), pred43_sbin4.at(i), D3_sbin4, A_sbin4, B_sbin4, weights,
                     lower.at(i), upper.at(i));
    }

    if(!plot_only){
      PrintLine(tex.at(i), Aval.at(i), Aup.at(i), Adown.at(i), Bval.at(i), Bup.at(i),
                Bdown.at(i), C3val.at(i), C3up.at(i), C3down.at(i), D3val.at(i),
                D3up.at(i), D3down.at(i), C2val.at(i), C2up.at(i), C2down.at(i),
                D2val.at(i), D2up.at(i), D2down.at(i), kappa3val.at(i), kappa3up.at(i),
                kappa3down.at(i), kappa2val.at(i), kappa2up.at(i), kappa2down.at(i),
                pred23.at(i), predup23.at(i), preddown23.at(i), pred24.at(i),
                predup24.at(i), preddown24.at(i), pred34.at(i), predup34.at(i),
                preddown34.at(i));
    }
  }

  if(plot_only){
    TCanvas canvas;
    TH1D h_closure_test("h_closure_test", "Closure Test;Sample;Events/19.4 fb^{-1}", 15, 0.5, 15.5);
    h_closure_test.GetXaxis()->SetBinLabel(1,"4b,all S");
    h_closure_test.GetXaxis()->SetBinLabel(2,"3b,all S");
    h_closure_test.GetXaxis()->SetBinLabel(3,"2b,all S");
    h_closure_test.GetXaxis()->SetBinLabel(4,"4b,S-bin 1");
    h_closure_test.GetXaxis()->SetBinLabel(5,"3b,S-bin 1");
    h_closure_test.GetXaxis()->SetBinLabel(6,"2b,S-bin 1");
    h_closure_test.GetXaxis()->SetBinLabel(7,"4b,S-bin 2");
    h_closure_test.GetXaxis()->SetBinLabel(8,"3b,S-bin 2");
    h_closure_test.GetXaxis()->SetBinLabel(9,"2b,S-bin 2");
    h_closure_test.GetXaxis()->SetBinLabel(10,"4b,S-bin 3");
    h_closure_test.GetXaxis()->SetBinLabel(11,"3b,S-bin 3");
    h_closure_test.GetXaxis()->SetBinLabel(12,"2b,S-bin 3");
    h_closure_test.GetXaxis()->SetBinLabel(13,"4b,S-bin 4");
    h_closure_test.GetXaxis()->SetBinLabel(14,"3b,S-bin 4");
    h_closure_test.GetXaxis()->SetBinLabel(15,"2b,S-bin 4");
    std::vector<double> val(0), up(0), down(0), halves(15,0.5), xval(0);
    xval.push_back(1.0); val.push_back(Aval.at(0)); up.push_back(Aup.at(0)); down.push_back(Adown.at(0));
    xval.push_back(2.0); val.push_back(C3val.at(0)); up.push_back(C3up.at(0)); down.push_back(C3down.at(0));
    xval.push_back(3.0); val.push_back(C2val.at(0)); up.push_back(C2up.at(0)); down.push_back(C2down.at(0));
    xval.push_back(4.0); val.push_back(Aval_sbin1.at(0)); up.push_back(Aup_sbin1.at(0)); down.push_back(Adown_sbin1.at(0));
    xval.push_back(5.0); val.push_back(C3val_sbin1.at(0)); up.push_back(C3up_sbin1.at(0)); down.push_back(C3down_sbin1.at(0));
    xval.push_back(6.0); val.push_back(C2val_sbin1.at(0)); up.push_back(C2up_sbin1.at(0)); down.push_back(C2down_sbin1.at(0));
    xval.push_back(7.0); val.push_back(Aval_sbin2.at(0)); up.push_back(Aup_sbin2.at(0)); down.push_back(Adown_sbin2.at(0));
    xval.push_back(8.0); val.push_back(C3val_sbin2.at(0)); up.push_back(C3up_sbin2.at(0)); down.push_back(C3down_sbin2.at(0));
    xval.push_back(9.0); val.push_back(C2val_sbin2.at(0)); up.push_back(C2up_sbin2.at(0)); down.push_back(C2down_sbin2.at(0));
    xval.push_back(10.0); val.push_back(Aval_sbin3.at(0)); up.push_back(Aup_sbin3.at(0)); down.push_back(Adown_sbin3.at(0));
    xval.push_back(11.0); val.push_back(C3val_sbin3.at(0)); up.push_back(C3up_sbin3.at(0)); down.push_back(C3down_sbin3.at(0));
    xval.push_back(12.0); val.push_back(C2val_sbin3.at(0)); up.push_back(C2up_sbin3.at(0)); down.push_back(C2down_sbin3.at(0));
    xval.push_back(13.0); val.push_back(Aval_sbin4.at(0)); up.push_back(Aup_sbin4.at(0)); down.push_back(Adown_sbin4.at(0));
    xval.push_back(14.0); val.push_back(C3val_sbin4.at(0)); up.push_back(C3up_sbin4.at(0)); down.push_back(C3down_sbin4.at(0));
    xval.push_back(15.0); val.push_back(C2val_sbin4.at(0)); up.push_back(C2up_sbin4.at(0)); down.push_back(C2down_sbin4.at(0));
    TGraphAsymmErrors g_closure_test(xval.size(), &xval.at(0), &val.at(0), &halves.at(0), &halves.at(0), &up.at(0), &down.at(0));  
    g_closure_test.SetFillColor(1);
    g_closure_test.SetLineColor(1);
    g_closure_test.SetFillStyle(3003);
    
    std::vector<double> from2(0), from3(0), from4(0);
    std::vector<double> from2_x(0), from3_x(0), from4_x(0);
    std::vector<double> from2_up(0), from3_up(0), from4_up(0);
    std::vector<double> from2_down(0), from3_down(0), from4_down(0);
    const double shift(0.125);
    from2_x.push_back(1.0-shift);
    from2_x.push_back(2.0-shift);
    from2_x.push_back(4.0-shift);
    from2_x.push_back(5.0-shift);
    from2_x.push_back(7.0-shift);
    from2_x.push_back(8.0-shift);
    from2_x.push_back(10.0-shift);
    from2_x.push_back(11.0-shift);
    from2_x.push_back(13.0-shift);
    from2_x.push_back(14.0-shift);
    from3_x.push_back(1.0+shift);
    from3_x.push_back(3.0-shift);
    from3_x.push_back(4.0+shift);
    from3_x.push_back(6.0-shift);
    from3_x.push_back(7.0+shift);
    from3_x.push_back(9.0-shift);
    from3_x.push_back(10.0+shift);
    from3_x.push_back(12.0-shift);
    from3_x.push_back(13.0+shift);
    from3_x.push_back(15.0-shift);
    from4_x.push_back(2.0+shift);
    from4_x.push_back(3.0+shift);
    from4_x.push_back(5.0+shift);
    from4_x.push_back(6.0+shift);
    from4_x.push_back(8.0+shift);
    from4_x.push_back(9.0+shift);
    from4_x.push_back(11.0+shift);
    from4_x.push_back(12.0+shift);
    from4_x.push_back(14.0+shift);
    from4_x.push_back(15.0+shift);
    from2.push_back(pred24.at(0));
    from2.push_back(pred23.at(0));
    from2.push_back(pred24_sbin1.at(0));
    from2.push_back(pred23_sbin1.at(0));
    from2.push_back(pred24_sbin2.at(0));
    from2.push_back(pred23_sbin2.at(0));
    from2.push_back(pred24_sbin3.at(0));
    from2.push_back(pred23_sbin3.at(0));
    from2.push_back(pred24_sbin4.at(0));
    from2.push_back(pred23_sbin4.at(0));
    from2_up.push_back(predup24.at(0));
    from2_up.push_back(predup23.at(0));
    from2_up.push_back(predup24_sbin1.at(0));
    from2_up.push_back(predup23_sbin1.at(0));
    from2_up.push_back(predup24_sbin2.at(0));
    from2_up.push_back(predup23_sbin2.at(0));
    from2_up.push_back(predup24_sbin3.at(0));
    from2_up.push_back(predup23_sbin3.at(0));
    from2_up.push_back(predup24_sbin4.at(0));
    from2_up.push_back(predup23_sbin4.at(0));
    from2_down.push_back(preddown24.at(0));
    from2_down.push_back(preddown23.at(0));
    from2_down.push_back(preddown24_sbin1.at(0));
    from2_down.push_back(preddown23_sbin1.at(0));
    from2_down.push_back(preddown24_sbin2.at(0));
    from2_down.push_back(preddown23_sbin2.at(0));
    from2_down.push_back(preddown24_sbin3.at(0));
    from2_down.push_back(preddown23_sbin3.at(0));
    from2_down.push_back(preddown24_sbin4.at(0));
    from2_down.push_back(preddown23_sbin4.at(0));
    from3.push_back(pred34.at(0));
    from3.push_back(pred32.at(0));
    from3.push_back(pred34_sbin1.at(0));
    from3.push_back(pred32_sbin1.at(0));
    from3.push_back(pred34_sbin2.at(0));
    from3.push_back(pred32_sbin2.at(0));
    from3.push_back(pred34_sbin3.at(0));
    from3.push_back(pred32_sbin3.at(0));
    from3.push_back(pred34_sbin4.at(0));
    from3.push_back(pred32_sbin4.at(0));
    from3_up.push_back(predup34.at(0));
    from3_up.push_back(predup32.at(0));
    from3_up.push_back(predup34_sbin1.at(0));
    from3_up.push_back(predup32_sbin1.at(0));
    from3_up.push_back(predup34_sbin2.at(0));
    from3_up.push_back(predup32_sbin2.at(0));
    from3_up.push_back(predup34_sbin3.at(0));
    from3_up.push_back(predup32_sbin3.at(0));
    from3_up.push_back(predup34_sbin4.at(0));
    from3_up.push_back(predup32_sbin4.at(0));
    from3_down.push_back(preddown34.at(0));
    from3_down.push_back(preddown32.at(0));
    from3_down.push_back(preddown34_sbin1.at(0));
    from3_down.push_back(preddown32_sbin1.at(0));
    from3_down.push_back(preddown34_sbin2.at(0));
    from3_down.push_back(preddown32_sbin2.at(0));
    from3_down.push_back(preddown34_sbin3.at(0));
    from3_down.push_back(preddown32_sbin3.at(0));
    from3_down.push_back(preddown34_sbin4.at(0));
    from3_down.push_back(preddown32_sbin4.at(0));
    from4.push_back(pred43.at(0));
    from4.push_back(pred42.at(0));
    from4.push_back(pred43_sbin1.at(0));
    from4.push_back(pred42_sbin1.at(0));
    from4.push_back(pred43_sbin2.at(0));
    from4.push_back(pred42_sbin2.at(0));
    from4.push_back(pred43_sbin3.at(0));
    from4.push_back(pred42_sbin3.at(0));
    from4.push_back(pred43_sbin4.at(0));
    from4.push_back(pred42_sbin4.at(0));
    from4_up.push_back(predup43.at(0));
    from4_up.push_back(predup42.at(0));
    from4_up.push_back(predup43_sbin1.at(0));
    from4_up.push_back(predup42_sbin1.at(0));
    from4_up.push_back(predup43_sbin2.at(0));
    from4_up.push_back(predup42_sbin2.at(0));
    from4_up.push_back(predup43_sbin3.at(0));
    from4_up.push_back(predup42_sbin3.at(0));
    from4_up.push_back(predup43_sbin4.at(0));
    from4_up.push_back(predup42_sbin4.at(0));
    from4_down.push_back(preddown43.at(0));
    from4_down.push_back(preddown42.at(0));
    from4_down.push_back(preddown43_sbin1.at(0));
    from4_down.push_back(preddown42_sbin1.at(0));
    from4_down.push_back(preddown43_sbin2.at(0));
    from4_down.push_back(preddown42_sbin2.at(0));
    from4_down.push_back(preddown43_sbin3.at(0));
    from4_down.push_back(preddown42_sbin3.at(0));
    from4_down.push_back(preddown43_sbin4.at(0));
    from4_down.push_back(preddown42_sbin4.at(0));
    h_closure_test.SetStats(0);
    h_closure_test.SetLineColor(0);
    
    std::vector<double> zeroes(from2.size(), 0.0);
    
    TGraphAsymmErrors graph2(from2_x.size(), &from2_x.at(0), &from2.at(0), &zeroes.at(0), &zeroes.at(0), &from2_up.at(0), &from2_down.at(0));
    TGraphAsymmErrors graph3(from3_x.size(), &from3_x.at(0), &from3.at(0), &zeroes.at(0), &zeroes.at(0), &from3_up.at(0), &from3_down.at(0));
    TGraphAsymmErrors graph4(from4_x.size(), &from4_x.at(0), &from4.at(0), &zeroes.at(0), &zeroes.at(0), &from4_up.at(0), &from4_down.at(0));
    graph2.SetLineColor(2);
    graph3.SetLineColor(3);
    graph4.SetLineColor(4);
    graph2.SetMarkerStyle(20);
    graph3.SetMarkerStyle(20);
    graph4.SetMarkerStyle(20);
    graph2.SetMarkerColor(2);
    graph3.SetMarkerColor(3);
    graph4.SetMarkerColor(4);
    
    double max(TMath::MaxElement(g_closure_test.GetN(), g_closure_test.GetY()));
    double max2(TMath::MaxElement(graph2.GetN(), graph2.GetY()));
    double max3(TMath::MaxElement(graph3.GetN(), graph3.GetY()));
    double max4(TMath::MaxElement(graph4.GetN(), graph4.GetY()));
    if(graph2.GetMaximum()>max) max=graph2.GetMaximum();
    if(graph3.GetMaximum()>max) max=graph3.GetMaximum();
    if(graph4.GetMaximum()>max) max=graph4.GetMaximum();
    if(max2>max) max=max2;
    if(max3>max) max=max3;
    if(max4>max) max=max4;
    h_closure_test.SetMaximum(1.1*max);
  
    h_closure_test.Draw();
    g_closure_test.Draw("2same");
    g_closure_test.Draw("psame");
    graph2.Draw("psame");
    graph3.Draw("psame");
    graph4.Draw("psame");
    
    TLegend legend(0.8,0.85,1.0,1.0);
    legend.AddEntry(&g_closure_test,"Observed","lpef");
    legend.AddEntry(&graph2,"Pred. from 2b","lpe");
    legend.AddEntry(&graph3,"Pred. from 3b","lpe");
    legend.AddEntry(&graph4,"Pred. from 4b","lpe");
    legend.Draw("same");

    if(mc_plot){
      canvas.Print("plots/mc_closure_test.pdf");
      canvas.SetLogy(1);
      canvas.Print("plots/mc_closure_test_log.pdf");
    }else{
      canvas.Print("plots/data_closure_test.pdf");
      canvas.SetLogy(1);
      canvas.Print("plots/data_closure_test_log.pdf");
    }
  }
  
  KillTrees(trees);
  KillFiles(files);
}

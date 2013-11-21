#include <iostream>
#include <string>
#include <unistd.h>
#include "TH1.h"
#include "event_handler.hpp"
#include "weights.hpp"
#include <ctime>

using namespace std;

int main(int argc, char *argv[]){
  time_t startTime, curTime;
  time(&startTime);

  TH1::SetDefaultSumw2(true);
  std::string inFilename("");
  bool iscfA(true), isfast(true);
  int c(0), Nentries(0);
  while((c=getopt(argc, argv, "n:i:c:f"))!=-1){
    switch(c){
    case 'n':
      Nentries=atoi(optarg);
      break;
    case 'i':
      inFilename=optarg;
      break;
    case 'c':
      iscfA=false;
      break;
    case 'f':
      isfast=false;
      break;
    }
  }
  
  std::string outFilename("");
  if(iscfA){
    outFilename="raw_plots_and_values/"+inFilename+".root";
    inFilename="/net/cms2/cms2r0/cfA/"+inFilename+"/cfA_"+inFilename+"*.root";
  }else{
    std::string baseName(inFilename);
    size_t pos(baseName.find(".root"));
    if(pos!=std::string::npos){
      baseName.erase(pos);
    }
    pos=baseName.rfind("/");
    if(pos!=std::string::npos){
      if(pos!=baseName.size()-1){
        baseName.erase(0,pos+1);
      }else{
        baseName.append("file_name_ended_with_slash");
      }
    }
    outFilename="raw_plots_and_values/"+baseName+".root";
    std::cout << inFilename << "\n" << baseName << "\n" << outFilename << "\n";
  }

  WeightCalculator w(19399);
  EventHandler eH(inFilename, false, w.GetWeight(inFilename), isfast); 

  time(&curTime);
  cout<<"Getting started takes "<<difftime(curTime,startTime)<<" seconds"<<endl;
  eH.MakePlots13Tev(outFilename, Nentries);

  time(&curTime);
  cout<<Nentries<<" events took "<<difftime(curTime,startTime)<<" seconds"<<endl;
}

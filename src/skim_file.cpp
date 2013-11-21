//Skims MC backgrounds using the full path to a local .root file (/path/to/local/file.root) or cfA sample name (TTJets_blahblahblah)

#include <string>
#include <iostream>
#include <unistd.h>
#include "event_handler.hpp"
#include "weights.hpp"

int main(int argc, char *argv[]){
  std::string inFilename("");
  bool isLocal(false);
  int c(0);
  int chargino_mass(-1), LSP_mass(-1);
  std::string chargino_mass_string(""), LSP_mass_string("");
  while((c=getopt(argc, argv, "i:M:m:l"))!=-1){
    switch(c){
    case 'i':
      inFilename=optarg;
      break;
    case 'M':
      chargino_mass_string=optarg;
      chargino_mass=atoi(optarg);
      break;
    case 'm':
      LSP_mass_string=optarg;
      LSP_mass=atoi(optarg);
      break;
    case 'l':
      isLocal=true;
      break;
    }
  }

  std::string outFilename("");
  if(!isLocal){
    outFilename="../data/"+inFilename+"_SyncSkim.root";
    inFilename="/net/cms2/cms2r0/cfA/"+inFilename+"/cfA_"+inFilename+"*.root";
  }else{
    std::string baseName(inFilename);
    size_t pos=baseName.find(".root");
    if(pos!=std::string::npos){
      baseName.erase(pos);
    }
    pos=baseName.rfind("/");
    if(pos!=std::string::npos){
      baseName.erase(0,pos+1);
    }
    outFilename="../data/"+baseName+"_SyncSkim.root";
  }

  if(chargino_mass>=0 || LSP_mass>=0){
    if(chargino_mass>=0 && LSP_mass>=0){
      std::string outFilename_tmp(outFilename);
      bool failed(false);
      const std::string::size_type chargino_pos(outFilename_tmp.find("mChargino"));
      const std::string::size_type chargino_start(outFilename_tmp.find("-",chargino_pos));
      const std::string::size_type chargino_end(outFilename_tmp.find("_",chargino_pos));
      if(chargino_pos!=std::string::npos
         && chargino_start!=std::string::npos
         && chargino_end!=std::string::npos){
        const std::string::size_type chargino_delta(chargino_end-chargino_start-1);
        outFilename_tmp.replace(chargino_start+1,chargino_delta,chargino_mass_string);
      }else{
        failed=true;
      }
      const std::string::size_type LSP_pos(outFilename_tmp.find("mLSP"));
      const std::string::size_type LSP_start(outFilename_tmp.find("-",LSP_pos));
      const std::string::size_type LSP_end(outFilename_tmp.find("_",LSP_pos));
      if(LSP_pos!=std::string::npos
         && LSP_start!=std::string::npos
         && LSP_end!=std::string::npos){
        const std::string::size_type LSP_delta(LSP_end-LSP_start-1);
        outFilename_tmp.replace(LSP_start+1,LSP_delta,LSP_mass_string);
      }else{
        failed=true;
      }
      if(!failed){
        outFilename=outFilename_tmp;
      }
    }else{
      std::cerr << "Error: Must specify both chargino and LSP mass or neither; cannot specify only one." << std::endl;
      chargino_mass=-1;
      LSP_mass=-1;
    }
  }

  WeightCalculator w(19399);
  EventHandler eH(inFilename, false, w.GetWeight(outFilename));
  eH.Skim(outFilename, chargino_mass, LSP_mass);
}

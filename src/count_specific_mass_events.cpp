#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <unistd.h>
#include "TFile.h"
#include "TChain.h"

int main(int argc, char *argv[]){
  std::string filename("");
  unsigned int small_mass(0), big_mass(0);
  char opt(' ');
  while(( opt=getopt(argc, argv, "f:m:M:") )!=-1){
    switch(opt){
    case 'f':
      filename=optarg;
      break;
    case 'm':
      small_mass=atoi(optarg);
      break;
    case 'M':
      big_mass=atoi(optarg);
      break;
    default:
      std::cerr << "Error in " << argv[0] << ": '" << opt
                << "' is not a valid option." << std::endl;
    }
  }

  TFile file(filename.c_str(),"read");
  if(file.IsOpen()){
    if(!file.IsZombie()){
      TChain *tree(static_cast<TChain*>(file.Get("configurableAnalysis/eventB")));
      if(tree!=NULL){
        std::string *model_params(NULL);
        tree->SetBranchStatus("*",0);
        tree->SetBranchStatus("model_params",1);
        tree->SetBranchAddress("model_params", &model_params);
        unsigned int counter(0);
        const unsigned int nentries(tree->GetEntries());
        for(unsigned int entry(0); entry<nentries; ++entry){
          tree->GetEntry(entry);
          if(false && small_mass && big_mass){}//Hardcode mass for now...need to fix
          if(model_params->find("chargino300")!=std::string::npos
             && model_params->find("bino1_")!=std::string::npos){
            std::cout << model_params->c_str() << std::endl;
            ++counter;
          }
        }
        std::cout << counter << std::endl;
      }
    }
    file.Close();
  }
}

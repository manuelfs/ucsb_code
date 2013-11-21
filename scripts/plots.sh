#! /bin/bash

export Nevents=10000

#./scripts/make_plots.exe -i SMS-T1tttt_2J_mGo-825_mLSP-225_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V19-v1_AODSIM_UCSB1788reshuf_v68 -n $Nevents
#./scripts/make_plots.exe -i SMS-T1tttt_2J_mGo-845to3000_mLSP-1to1355_TuneZ2star_14TeV-madgraph-tauola_Summer12-START53_V7C_FSIM_PU_S12-v1_AODSIM_UCSB1949reshuf_v71  -n $Nevents

./scripts/make_plots.exe -i TTbar_TuneZ2star_13TeV-pythia6-tauola_Summer13dr53X-PU45bx25_START53_V19D-v2_AODSIM_UCSB1933_v71 -n $Nevents
./scripts/make_plots.exe -i TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71 -n $Nevents
./scripts/make_plots.exe -i TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71 -n $Nevents
./scripts/make_plots.exe -i TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71 -n $Nevents







wait

echo "Done"
exit 0;
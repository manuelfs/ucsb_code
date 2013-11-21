#! /bin/bash

./scripts/skim_file.exe -i SMS-T1tttt_2J_mGo-825_mLSP-225_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V19-v1_AODSIM_UCSB1788reshuf_v68 &> logs/make_8TeV-T1tttt.log &
./scripts/skim_file.exe -i SMS-T1tttt_2J_mGo-845to3000_mLSP-1to1355_TuneZ2star_14TeV-madgraph-tauola_Summer12-START53_V7C_FSIM_PU_S12-v1_AODSIM_UCSB1949reshuf_v71 &> logs/make_13TeV-T1tttt.log &

./scripts/skim_file.exe -i TTbar_TuneZ2star_13TeV-pythia6-tauola_Summer13dr53X-PU45bx25_START53_V19D-v2_AODSIM_UCSB1933_v71 &> logs/make_13TeV_TTbar.log &
./scripts/skim_file.exe -i TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71 &> logs/make_skim_1883.log &
./scripts/skim_file.exe -i TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71 &> logs/make_skim_1880.log &
./scripts/skim_file.exe -i TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71 &> logs/make_skim_1884.log &







wait

echo "Done"
exit 0;

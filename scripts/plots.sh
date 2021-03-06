#! /bin/bash

export Nevents=1000000

# Signal @14TeV: 40,321 entries
./scripts/make_plots.exe -i SMS-T1tttt_2J_mGo-845to3000_mLSP-1to1355_TuneZ2star_14TeV-madgraph-tauola_Summer12-START53_V7C_FSIM_PU_S12-v1_AODSIM_UCSB1949reshuf_v71  -n $Nevents -m 1145_800

# Signal @14TeV:
./scripts/make_plots.exe -i SMS-T1tttt_2J_mGo-845to3000_mLSP-1to1355_TuneZ2star_14TeV-madgraph-tauola_Summer12-START53_V7C_FSIM_PU_S12-v1_AODSIM_UCSB1949reshuf_v71  -n $Nevents -m 1145_500
# Signal @8TeV: 127,823 entries
./scripts/make_plots.exe -i  SMS-MadGraph_Pythia6Zstar_8TeV_T1tttt_2J_mGo-1100to1400_mLSP-525to1000_25GeVX25GeV_Binning_Summer12-START52_V9_FSIM-v2_AODSIM_UCSB1739reshuf_v68   -n $Nevents -m 1150_800_

# All ttbar @13TeV: 997,120 entries
./scripts/make_plots.exe -i TTbar_TuneZ2star_13TeV-pythia6-tauola_Summer13dr53X-PU45bx25_START53_V19D-v2_AODSIM_UCSB1933_v71 -n $Nevents
# Leptonic ttbar @13TeV: 12,011,428 entries
./scripts/make_plots.exe -i TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71 -n $Nevents
# Hadronic ttbar @13TeV: 31,223,821 entries
./scripts/make_plots.exe -i TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71 -n $Nevents
# Semi-leptonic ttbar @13TeV: 24,953,451 entries
./scripts/make_plots.exe -i TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71 -n $Nevents


wait

echo "Done"
exit 0;

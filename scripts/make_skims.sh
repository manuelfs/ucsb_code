#! /bin/bash

./scripts/skim_file.exe -i WbbJetsToLNu_Massive_TuneZ2star_8TeV-madgraph-pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1859_v71 &> logs/make_skim_1859.log &

./scripts/skim_file.exe -i TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71 &> logs/make_skim_1883.log &
./scripts/skim_file.exe -i TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71 &> logs/make_skim_1880.log &
./scripts/skim_file.exe -i TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71 &> logs/make_skim_1884.log &

./scripts/skim_file.exe -i TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1857_v71 &> logs/make_skim_1857.log &
./scripts/skim_file.exe -i TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1856_v71 &> logs/make_skim_1856.log &

./scripts/skim_file.exe -i ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1887_v71 &> logs/make_skim_1887.log &
./scripts/skim_file.exe -i ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1889_v71 &> logs/make_skim_1889.log &
./scripts/skim_file.exe -i ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1888_v71 &> logs/make_skim_1888.log &
./scripts/skim_file.exe -i ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1891_v71 &> logs/make_skim_1891.log &
./scripts/skim_file.exe -i ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1890_v71 &> logs/make_skim_1890.log &

./scripts/skim_file.exe -i ZH_ZToBB_HToBB_M-125_8TeV-powheg-herwigpp_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1868_v71 &> logs/make_skim_1868.log &
./scripts/skim_file.exe -i WH_WToLNu_HToBB_M-125_8TeV-powheg-herwigpp_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1858_v71 &> logs/make_skim_1858.log &

./scripts/skim_file.exe -i WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1874_v71 &> logs/make_skim_1874.log &
./scripts/skim_file.exe -i ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1876_v71 &> logs/make_skim_1876.log &

./scripts/skim_file.exe -i TTH_HToBB_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1855_v71 &> logs/make_skim_1855.log &

./scripts/skim_file.exe -i QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM_UCSB1897_v71 &> logs/make_skim_1897.log &
./scripts/skim_file.exe -i QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1898_v71 &> logs/make_skim_1898.log &
./scripts/skim_file.exe -i QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1899_v71 &> logs/make_skim_1899.log &
./scripts/skim_file.exe -i QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1900_v71 &> logs/make_skim_1900.log &
./scripts/skim_file.exe -i QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1901_v71 &> logs/make_skim_1901.log &
./scripts/skim_file.exe -i QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1902_v71 &> logs/make_skim_1902.log &
./scripts/skim_file.exe -i QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1903_v71 &> logs/make_skim_1903.log &
./scripts/skim_file.exe -i QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1904_v71 &> logs/make_skim_1904.log &
./scripts/skim_file.exe -i QCD_Pt-1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1905_v71 &> logs/make_skim_1905.log &

./scripts/skim_file.exe -i BJets_HT-250To500_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1893_v71 &> logs/make_skim_1893.log &
./scripts/skim_file.exe -i BJets_HT-500To1000_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1894_v71 &> logs/make_skim_1894.log &
./scripts/skim_file.exe -i BJets_HT-1000ToInf_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1895_v71 &> logs/make_skim_1895.log &

./scripts/skim_file.exe -i W2JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1877_v71 &> logs/make_skim_1877.log &
./scripts/skim_file.exe -i W3JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1878_v71 &> logs/make_skim_1878.log &
./scripts/skim_file.exe -i W4JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1879_v71 &> logs/make_skim_1879.log &

./scripts/skim_file.exe -i T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1860_v71 &> logs/make_skim_1860.log &
./scripts/skim_file.exe -i T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1861_v71 &> logs/make_skim_1861.log &
./scripts/skim_file.exe -i T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1862_v71 &> logs/make_skim_1862.log &
./scripts/skim_file.exe -i Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1864_v71 &> logs/make_skim_1864.log &
./scripts/skim_file.exe -i Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1865_v71 &> logs/make_skim_1865.log &
./scripts/skim_file.exe -i Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1866_v71 &> logs/make_skim_1866.log &

./scripts/skim_file.exe -i MET_Run2012A-13Jul2012-v1_AOD_UCSB1852_v71 &> logs/make_skim_1852.log &
./scripts/skim_file.exe -i MET_Run2012B-13Jul2012-v1_AOD_UCSB1853_v71 &> logs/make_skim_1853.log &
./scripts/skim_file.exe -i MET_Run2012C-24Aug2012-v1_AOD_UCSB1854_v71 &> logs/make_skim_1854.log &
./scripts/skim_file.exe -i MET_Run2012C-PromptReco-v2_AOD_UCSB1867_v71 &> logs/make_skim_1867.log &
./scripts/skim_file.exe -i MET_Run2012D-PromptReco-v1_AOD_UCSB1869_v71 &> logs/make_skim_1869.log &
./scripts/skim_file.exe -i MET_Run2012D-PromptReco-v1_AOD_UCSB1870_v71 &> logs/make_skim_1870.log &

./scripts/skim_file.exe -i SMS-TChiHH_2b2b_2J_mChargino-130to325_mLSP-1to195_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872_v71 -M 200 -m 1 &> logs/make_skim_1872_200_1.log &
./scripts/skim_file.exe -i SMS-TChiHH_2b2b_2J_mChargino-130to325_mLSP-1to195_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872_v71 -M 250 -m 1 &> logs/make_skim_1872_250_1.log &
./scripts/skim_file.exe -i SMS-TChiHH_2b2b_2J_mChargino-130to325_mLSP-1to195_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872_v71 -M 300 -m 1 &> logs/make_skim_1872_300_1.log &
./scripts/skim_file.exe -i SMS-TChiHH_2b2b_2J_mChargino-350to500_mLSP-1to370_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1871_v71 -M 350 -m 1 &> logs/make_skim_1871_350_1.log &
./scripts/skim_file.exe -i SMS-TChiHH_2b2b_2J_mChargino-350to500_mLSP-1to370_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1871_v71 -M 400 -m 1 &> logs/make_skim_1871_400_1.log &
./scripts/skim_file.exe -i SMS-TChiHH_2b2b_2J_mChargino-350to500_mLSP-1to370_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1871_v71 -M 450 -m 1 &> logs/make_skim_1871_450_1.log &

wait

echo "Done"
exit 0;

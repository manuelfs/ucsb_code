#! /bin/bash

cp src/event_handler.cpp testing.cpp

echo "Running without top pt reweighting."
sed 's/true \&\& isttbar?GetTopPtWeight():1.0/false \&\& isttbar?GetTopPtWeight():1.0/g' -i src/event_handler.cpp
make clean
./compile.sh
./scripts/make_plots.sh &> logs/make_plots.log
rm -rf raw_plots_and_values_no_pt_weight
cp -r raw_plots_and_values raw_plots_and_values_no_pt_weight
./scripts/calc_abcd.exe &> logs/calc_abcd_no_ptweight.txt
sed 's/false \&\& isttbar?GetTopPtWeight():1.0/true \&\& isttbar?GetTopPtWeight():1.0/g' -i src/event_handler.cpp

echo "Running with increased gluon splitting."
sed 's/HasGluonSplitting()?1.0:1.0/HasGluonSplitting()?1.5:1.0/g' -i src/event_handler.cpp
make clean
./compile.sh
./scripts/make_plots.sh &> logs/make_plots.log
rm -rf raw_plots_and_values_gluon_up
cp -r raw_plots_and_values raw_plots_and_values_gluon_up
./scripts/calc_abcd.exe &> logs/calc_abcd_gluon_up.txt
sed 's/HasGluonSplitting()?1.5:1.0/HasGluonSplitting()?1.0:1.0/g' -i src/event_handler.cpp

echo "Running with decreased gluon splitting."
sed 's/HasGluonSplitting()?1.0:1.0/HasGluonSplitting()?0.5:1.0/g' -i src/event_handler.cpp
make clean
./compile.sh
./scripts/make_plots.sh &> logs/make_plots.log
rm -rf raw_plots_and_values_gluon_down
cp -r raw_plots_and_values raw_plots_and_values_gluon_down
./scripts/calc_abcd.exe &> logs/calc_abcd_gluon_down.txt
sed 's/HasGluonSplitting()?0.5:1.0/HasGluonSplitting()?1.0:1.0/g' -i src/event_handler.cpp

echo "Running normally."
make clean
./compile.sh
./scripts/make_plots.sh &> logs/make_plots.log
./scripts/calc_abcd.exe &> logs/calc_abcd.log

echo "diff test:"
diff src/event_handler.cpp testing.cpp
rm -rf testing.cpp
echo "Done."


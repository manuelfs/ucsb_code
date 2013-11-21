#! /bin/bash

shopt -s nullglob

for file in ../data/*.root
do
  noext=${file%.root}
  ./scripts/make_plots.exe -i "$file" &> "logs/make_plots_${noext##*/}.log" &
done

wait
echo "Done"
exit 0;

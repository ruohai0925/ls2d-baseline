#!/bin/bash
# usage: run_table.sh <base inputs> <outprefix> "<extra key=val ...>" grid1 grid2 ...
base=$1; pre=$2; extra=$3; shift 3
S=out/_inputs; mkdir -p $S
for n in "$@"; do
  f=$S/inputs.$pre.$n; grep -v -E '^amr.n_cell|^out_dir' $base > $f; echo "amr.n_cell = $n $n" >> $f; echo "out_dir = out/$pre$n" >> $f
  for kv in $extra; do k=${kv%%=*}; grep -v "^$k " $f > $f.t; mv $f.t $f; echo "$k = ${kv#*=}" >> $f; done
  ./ls2d $f > out/log.$pre$n 2>&1 &
done
wait
printf "| grid | area | %% loss | err/L | grad_err_max |\n|---|---|---|---|---|\n"
for n in "$@"; do s=out/$pre$n/summary.txt; printf "| %s | %.4g | %.3g | %.3g | %.3g |\n" $n $(grep area_sharp_final $s|cut -d' ' -f2) $(grep area_loss_pct $s|cut -d' ' -f2) $(grep over_L $s|cut -d' ' -f2) $(grep grad_err_max $s|cut -d' ' -f2); done

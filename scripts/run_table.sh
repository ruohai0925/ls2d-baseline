#!/bin/bash
# Grid sweep for a flow input file: runs every grid in parallel and prints the convergence table.
# usage: run_table.sh <base inputs> <outprefix> "<extra key=val ...>" grid1 grid2 ...
#   e.g. scripts/run_table.sh tests/inputs.taylorgreen tg "" 32 64 128
base=$1; pre=$2; extra=$3; shift 3
S=out/_inputs; mkdir -p $S
for n in "$@"; do
  f=$S/inputs.$pre.$n; grep -v -E '^amr.n_cell|^out_dir' $base > $f; echo "amr.n_cell = $n $n" >> $f; echo "out_dir = out/$pre$n" >> $f
  for kv in $extra; do k=${kv%%=*}; grep -v "^$k " $f > $f.t; mv $f.t $f; echo "$k = ${kv#*=}" >> $f; done
  ./ls2d $f > out/log.$pre$n 2>&1 &
done
wait
printf "| grid | L2(u) error | max nodal divergence | steps | wall time (s) |\n|---|---|---|---|---|\n"
for n in "$@"; do s=out/$pre$n/summary.txt
  printf "| %s | %s | %s | %s | %s |\n" $n "$(grep '^L2_error_u' $s|cut -d' ' -f2)" "$(grep '^max_div_nodal' $s|cut -d' ' -f2)" \
    "$(grep '^steps' $s|cut -d' ' -f2)" "$(grep '^wall_time_s' $s|cut -d' ' -f2)"
done

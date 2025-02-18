#!/bin/bash

for k in {1..1000}
do
  for i in {1..18}
  do sleep 5
    Rscript estsim_ursus.R &
  done
  wait
done
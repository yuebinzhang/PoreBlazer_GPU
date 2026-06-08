#!/bin/bash
##$ -R y
#$ -l h_rt=47:00:00
#$ -l h_vmem=24G
#$ -j y
#$ -m a
#$ -cwd

#/bin/date
. /etc/profile
port=$((JOB_ID % 5000 + 20000))

ulimit -v
time ../poreblazer.exe < input.dat > results


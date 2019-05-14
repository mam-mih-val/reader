#!/bin/bash

hadesroot=/cvmfs/hades.gsi.de/install/6.12.06/hydra2-4.9w/defall.sh 

current_dir = $(pwd)
partition = main
time=8:00:00
executable=../build/DT_Reader

file_list =     ${1}
cmd =           ${2}
config =        ${3}
output_dir =    ${4}

log_dir=${output_dir}/log

mkdir -p $output_dir
mkdir -p $log_dir

n_runs=$(cat file_list | wc -l)

job_range=1-$n_runs

executable=${current_dir}/../build/DT_Reader

echo executable=$executable
echo output_dir=$output_dir
echo log_dir=$log_dir
echo n_runs=$n_runs
echo job_range=$job_range

sbatch -J DT_Reader -p $partition -t $time -a $job_range -e ${log_dir}/%A_%a.e -o ${log_dir}/%A_%a.o --export=executable=$executable,output_dir=$output_dir,file_list=$file_list,hadesroot=$hadesroot,cmd=$cmd,config=$config batch_run.sh
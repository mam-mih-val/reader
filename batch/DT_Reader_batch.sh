#!/bin/bash

hadesroot=/cvmfs/hades.gsi.de/install/6.12.06/hydra2-4.9w/defall.sh 

current_dir=$(pwd)
partition=main
time=8:00:00
executable=../build/DT_Reader
signal=adc
channelSelection=0
minSignal=0
maxSignal=999
pidCode=14

while [ "$#" -gt "2" ]
do
    #echo $#
    case "$1" in
    --flow) cmd=flow 
    echo found cmd=$cmd ;;
    --eventqa) cmd=eventqa
    echo found cmd=$cmd ;;
    --trackqa) cmd=trackqa
    echo found cmd=$cmd ;;
    --qvector) cmd=qvector
    echo found cmd=$cmd ;;
    --signal) signal=$2
    echo found signal=$signal
    shift ;;
    --perchannel) channelSelection=$2
    echo found Channel Selection = $channelSelection
    shift ;;
    --min) minSignal=$2
    echo found minimal signal = $minSignal
    shift ;;
    --max) maxSignal=$2
    echo found maximal signal = $maxSignal
    shift ;;
    --pid) pidCode=$2
    echo found pid code = $pidCode
    --) shift
    break ;;
    # *) echo Error: incorrect option
    esac
    shift
done
file_list=${1}
output_dir=${2}

log_dir=${output_dir}/log

mkdir -p $output_dir
mkdir -p $log_dir

n_runs=$(cat $file_list | wc -l)

job_range=1-$n_runs

executable=${current_dir}/../build/DT_Reader

echo file list=$file_list
echo executable=$executable
echo output_dir=$output_dir
echo log_dir=$log_dir
echo n_runs=$n_runs
echo job_range=$job_range
echo
echo Configuration:
echo cmd=$cmd
echo signal=$signal
echo Channel Selection=$channelSelection
echo Minimal Signal=$minSignal
echo Maximal Signal=$maxSignal
echo pid code=$pidCode

sbatch -J DT_Reader -p $partition -t $time -a $job_range -e ${log_dir}/%A_%a.e -o ${log_dir}/%A_%a.o --export=executable=$executable,output_dir=$output_dir,file_list=$file_list,hadesroot=$hadesroot,cmd=$cmd,signal=$signal,channelSelection=$channelSelection,minSignal=$minSignal,maxSignal=$maxSignal,pidCode=$pidCode batch_run.sh
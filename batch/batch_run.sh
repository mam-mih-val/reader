#!/bin/bash

format='+%Y/%m/%d-%H:%M:%S'

date $format

job_num=$(($SLURM_ARRAY_TASK_ID))
output_file=$output_dir/$job_num.root
input_file=`sed -n "${job_num}p" < $file_list`

echo "loading " $hadesroot
source $hadesroot
echo "executing $executable  $input_files $output_file $n_events"

$executable $cmd --signal $signal --perchannel $perchannel --min $minSignal --max $maxSignal $input_file $output_file

echo JOB FINISHED!
date $format
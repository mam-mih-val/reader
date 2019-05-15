#!/bin/bash

format='+%Y/%m/%d-%H:%M:%S'

date $format

job_num=$(($SLURM_ARRAY_TASK_ID))
output_file=$output_dir/$job_num.root
input_file=`sed "${job_num}q;d" $file_list`

echo "loading " $hadesroot
source $hadesroot
echo "executing $executable  $input_files $output_file $n_events"
$executable  $input_file $cmd $config $output_file

echo JOB FINISHED!
date $format
#!/bin/bash

format='+%Y/%m/%d-%H:%M:%S'

date $format

job_num=$(($SLURM_ARRAY_TASK_ID))
output_file=$output_dir/$job_num.root
input_file=`sed -n "${job_num}p" < $file_list`

echo "loading " $hadesroot
source $hadesroot

echo "executing $executable $cmd --signal $signal --perchannel $channelSelection --min $minSignal --max $maxSignal $input_file $output_file"

$executable $cmd --signal $signal --perchannel $channelSelection --min $minSignal --max $maxSignal $input_file $output_file

echo JOB FINISHED!
date $format
#!/bin/bash

start_time=$(date -R) 

#pipFold="$1"

pipFold="../TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER"

echo ">>> START CLEANING: $pipFold"

#all_datasets=( $(ls -d $pipFold/*) )

all_datasets=( 
"../TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER/TCGAcrc_msi_mss"
"../TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER/TCGAacc_acc_mutCTNNB1"
)

for datasetFold in "${all_datasets[@]}"; do
    echo $datasetFold


    all_log_files=( $(ls -t $datasetFold/*logFile.txt) )

    #seen_logFiles=()

    echo "NEED TO CHECK HOW UNDECLARE A VARIABLE"   
    #exit 
    unset seen_logFiles # unset will not work on declare read only try set instead of declare
    set -a seen_logFiles
    #declare -a seen_logFiles


    for logFile in "${all_log_files[@]}"; do
#        echo $logFile
        filename=`basename $logFile`
        echo $filename

        name_parts=( $(echo $filename | tr "_" "\n") )


        # if the first part of the name is not only number -> delete the file

        part1=${name_parts[0]}
        part2=${name_parts[1]}

        echo "name part1: $part1"
        echo "name part2: $part2"

        regex="^[0-9]+$"
        if [[ $part1 =~ $regex ]]; then
            # then look if part2 already seen
            # look if I already a file for this step -> if yes, delete the file, if no add to the array
            echo "find logFile for step: $part2"

            if [[ " ${seen_logFiles[@]} " =~ " ${part2} " ]]; then
                echo "this step is already in the array -> file should be deleted"
                echo rm -f $logFile
                rm -f $logFile
            else 
                echo "this step not yet in the array -> add it, do not delete the file"
                seen_logFiles=("${seen_logFiles[@]}" $part2)
            fi
        
        else
            "part1 does not match -> file should be deleted"
            echo rm -f $logFile
            rm -f $logFile
        fi

    done

#exit
done





###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

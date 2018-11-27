#!/bin/bash


caller="TopDom"
    
pipFold="/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_${caller}"

all_folders=( $(ls $pipFold/OUTPUT_FOLDER ))

for folder in ${all_folders[@]}; do 
    echo ">>> START: $folder"

    fpkm_file="$pipFold/OUTPUT_FOLDER/$folder/0_prepGeneData/rna_fpkmDT.Rdata"

    if [ ! -f $fpkm_file ] ; then
        echo "... $fpkm_file does not exist"
    fi


    permut_file="$pipFold/OUTPUT_FOLDER/$folder/5_runPermutationsMedian/permutationsDT.Rdata"

    if [ ! -f $permut_file ] ; then
        echo "... $permut_file does not exist"
    fi


    combined_file="$pipFold/OUTPUT_FOLDER/$folder/13_plotTopWilcoxTopCombined/venn_combined_wilcox_intersect_top50.svg"

    if [ ! -f $combined_file ] ; then
        echo "... $combined_file does not exist"
    fi


    auc_file="$pipFold/OUTPUT_FOLDER/$folder/170_score_auc_pval_withShuffle/allratio_auc_pval.Rdata"

    if [ ! -f $auc_file ] ; then
        echo "... $auc_file does not exist"
    fi


    nratios=`ls $pipFold/OUTPUT_FOLDER/$folder/8c_*/*.Rdata | wc -l `

    if [ $nratios -ne "6" ]; then
        echo "... missing files in 8c_ folder"
    fi




done

#!/bin/bash

pipfolder="/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2"
dataset="TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1"

all_datas=( "CL" "human" "EPZ" "GSK")

all_steps=( "14f2" "14i2")

picfile="allRatios_cumsum_obs_permut.svg"
convertFolder="tmp_pdf_pic"

for data in "${all_datas[@]}"; do
    outfolder="OUTPUT_FOLDER_EZH2_MAPQ_v2_$data"


    for step in "${all_steps[@]}"; do

        if [[ $step == "14i2" ]]; then 
            stepfolder="14i2_cumulAllDown_limited_AUC_randomTADsShuffle"
            patt="random_TADs"
        elif [[ $step == "14f2" ]]; then 
            stepfolder="14f2_cumulAllDown_limited_AUC"
            patt="random_genes"
        else
            echo "error"
            exit 1
        fi
#        mkdir -p $pipfolder/$outfolder/$convertFolder
#        echo inkscape -f $pipfolder/$outfolder/$dataset/$stepfolder/$picfile -A $pipfolder/$outfolder/$convertFolder/${data}_${step}_${patt}_cumsum.pdf
#        inkscape -f $pipfolder/$outfolder/$dataset/$stepfolder/$picfile -A $pipfolder/$outfolder/$convertFolder/${data}_${step}_${patt}_cumsum.pdf

        mkdir -p $pipfolder/$convertFolder
        echo inkscape -f $pipfolder/$outfolder/$dataset/$stepfolder/$picfile -A $pipfolder/$convertFolder/${data}_${step}_${patt}_cumsum.pdf
        inkscape -f $pipfolder/$outfolder/$dataset/$stepfolder/$picfile -A $pipfolder/$convertFolder/${data}_${step}_${patt}_cumsum.pdf

    done
done

exit 0















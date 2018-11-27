#!/usr/bin/bash

# running like
#  ./zzz_run_given_step_given_data.sh run_settings_TCGAcrc_entrez.R 14
#  ./zzz_run_given_step_given_data.sh <run_file> <step#>

scriptFolder="."

# ./zzz_run_given_step_given_data.sh ../TAD_DE_pipeline/run_settings_GSE40419_normal_cancer.R 0

set -e

args=( "$@" )
args_len=${#args[@]}

settingF=${args[0]}

i=1
while [[ $i -lt args_len ]]; do
    scriptF="$(ls $scriptFolder | grep -P ^${args[$i]}_.+R$)"
    cmd="Rscript $scriptFolder/$scriptF $settingF"
    echo "> $cmd"
    $cmd
    ((i++))
done


# ./zzz_run_given_step_given_data.sh BUILDDT_SETTING_FILES/run_settings_buildDT.R


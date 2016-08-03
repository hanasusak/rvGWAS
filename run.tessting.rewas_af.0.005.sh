#!/bin/sh

#$ -e /no_backup/GD/projects/CLL/germ_analysis/scripts/script_out
#$ -o /no_backup/GD/projects/CLL/germ_analysis/scripts/script_out

#$ -l h_rt=240:00:00 # time requested
#$ -l virtual_free=120G # RAM

#$ -q long-sl65
#$ -pe smp 5

#$ -m abe
#$ -M hana.susak@crg.eu

centR=/nfs/no_backup/GD/projects/CLL/germ_analysis/scripts/r/germline_risk_tests_v2.R

export R_MAX_MC_CORES=5
echo $R_MAX_MC_CORES

export READ_FROM='bash'
echo $READ_FROM

export R_SEED='160185'
echo $INPUT_ETA_START

export R_MIST='T'
echo $R_MIST
export R_KBAC='T'
echo $R_KBAC
export R_SKATO='T'
echo $R_SKATO

export R_COV_MAT_COLS=2
echo $R_COV_MAT_COLS

export R_PROJECT_COL='DB_Project'
echo $R_PROJECT_COL

export R_CASES_NAME='CLL'
echo $R_CASES_NAME

export R_PERMUTATIONS=100
echo $R_PERMUTATIONS

export R_CADD_TRH=10
echo $R_CADD_TRH

export R_CONTROLS_NUM=437
echo $R_CONTROLS_NUM

export R_PERM_ANSW='y'
echo $R_PERM_ANSW

MAT='/no_backup/GD/projects/CLL/germ_analysis/data/QC_output/mutations_filtered_samples_fix.txt'
samp_f='/no_backup/GD/projects/CLL/germ_analysis/data/QC_output/samples_info_pca.txt'
out='/no_backup/GD/projects/CLL/germ_analysis/temp_hana/mist_skato_kbac_out2'

Rscript $centR -m $MAT -d $samp_f -f 0.005 -o $out

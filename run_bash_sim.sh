#!/bin/sh

#$ -e /no_backup/GD/projects/RVAS/real_case_scenario/Simulation_2017-07-06-172351/Results_bash
#$ -o /no_backup/GD/projects/RVAS/real_case_scenario/Simulation_2017-07-06-172351/Results_bash

#$ -l h_rt=840:00:00 # time requested
#$ -l virtual_free=20G # RAM

#$ -q so-el7,long-sl7
#$ -pe smp 1

#$ -N Real_case_sim_RVAS


#$ -m abe
#$ -M hana.susak@crg.eu


## VARIABLES
centR=/nfs/no_backup/GD/projects/RVAS/scripts/RVAS_tests_prec.R

## basic params
export READ_FROM='bash'
echo $READ_FROM

export R_MAX_MC_CORES=1
echo $R_MAX_MC_CORES

export R_SEED='160185'
echo $R_SEED

export R_PERMUTATIONS=1
echo $R_PERMUTATIONS

## tests to run
export R_INLA='T'
echo $R_INLA
export R_MIST='T'
echo $R_MIST
export R_KBAC='T'
echo $R_KBAC
export R_SKATO='T'
echo $R_SKATO
export R_BURDEN='T'
echo $R_BURDEN
export MAF_LOGICAL='F'
echo $MAF_LOGICAL
export H_STEP= 0.001
echo $H_STEP
export R_BATI_SAVE='F'
echo $R_BATI_SAVE

## column names in samples file to be used as covariates
export R_COV_MAT_COLS='PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9'
echo $R_COV_MAT_COLS

## Column name in mutations file to be used as aggregation column
export R_AGG_COL='#Gene'
echo $R_AGG_COL

## Column name in samples file to be used as indicator for cases and columns 
export R_PROJECT_COL='case_ind'
echo $R_PROJECT_COL

export R_CASES_NAME='case'
echo $R_CASES_NAME

export R_CONTROLS_NAME='control'
echo $R_CONTROLS_NAME

## filtering params
export R_CHAR_FILT_COL='#ExonicFunction'
echo $R_CHAR_FILT_COL

export R_CHAR_FILT_VAL='synonymous SNV'
echo $R_CHAR_FILT_VAL

export R_NUM_FILT_COL='#EurEVSFrequency|#Eur1000GenomesFrequency|#ExAC_NFE|#Cadd2'
echo $R_NUM_FILT_COL

export R_NUM_FILT_VAL='h|h|h|l'
echo $R_NUM_FILT_VAL

export R_NUM_FILT_TRH='0.01|0.01|0.01|10'
echo $R_NUM_FILT_TRH

## setting random part variables
export R_DUMMY_VAR='#ExonicFunction'
echo $R_DUMMY_VAR

export R_DUMMY_MERGE='frameshift deletion, frameshift insertion, nonframeshift deletion, nonframeshift insertion | stopgain SNV, stoploss SNV, splicing'
echo $R_DUMMY_MERGE

export R_DUMMY_RENAME='indels|stop_gain_loss_splicing'
echo $R_DUMMY_RENAME

export R_NUMERIC_VAR='#Cadd2'
echo $R_NUMERIC_VAR

## perm settings
export R_CONTROLS_NUM=778
echo $R_CONTROLS_NUM

export R_PERM_ANSW='y'
echo $R_PERM_ANSW


MAT='/no_backup/GD/projects/RVAS/real_case_scenario/Simulation_2017-07-06-172351/sim.NO.muts.txt'
samp_f='/no_backup/GD/projects/RVAS/real_case_scenario/Simulation_2017-07-06-172351/sim.sample.info.txt'
out='/no_backup/GD/projects/RVAS/real_case_scenario/Simulation_2017-07-06-172351/Results_bash_no_sim_mut'

mkdir -p $out


Rscript $centR -m $MAT -d $samp_f -f 0.01 -o $out

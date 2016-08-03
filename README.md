# REWAS

---

In this repository are sequence of R scripts to perform REWAS (Rare-variant Exome Wide Association Study). To perform this analysis you will need R and installed nessesary packages (from R console instructions):
```R
install.packages(c("argparse", "ggplot2", "mixOmics","reshape2", "data.table", "MiST", "SKAT", "KBAC", "parallel"))
```

### Quality Control analysis
To start you need to perform inital QC check. For that you need at least two input files, multi-sample VCF-like file (with additional annotation, done by eDiVa) and Samples info file (tab separated, where Samples IDs are matching the one in multi-sample file). Additionaly you can add separete mutlis-sample file for InDels and samples IDs correction file (with mappings) in case some sample IDs are not matching in two before mentioned files.

For more information how to run this QC script type at your terminal:
```Shell
Rscript QC_and_filtering_multisaple_call_v2.R --h
```

### HUGO gene names correction
We observed that genes often have inconsistent names (e.g. sometimes in annotation is used old name or synonimous). We created script to correct this to the latest update of HUGO approved names, when possible, as there are ambiguous cases (e.g. one old name lead to 2 approved current names) which requires manual check up. Therefore we recomand to run this script to minimize possible confusion or errors with names.

For more information how to run this Correct names script type at your terminal:
```Shell
Rscript CorrectGeneNames.R --h
```
As database has restricted access, we recomand to just provide following:
* file path to multi-sample VCF-like file which was output from QC script
* column number in which HUGO gene names are in this file 
* optionally and recommandable for the first time you can print out steps of the script by setting -v T 

This will cause script to get mappings from url link, which is bit slower but is only nessesary to do once.


### Performing MiST, SKAT-O and KBAC tests
Now you are ready to perform tests. Your input files are output files from previous steps:
* VCF-like file with correcte HUGO names (output from Correct names script)
* Samples info file with PCA components added (output from QC script)
And additional parameters:
* Allele Frequncy threshold 
* List of genes you want to test (optional, if not set all genes in your seqence space will be used)
* Output folder (optional, if not set current directory will be used)

Other parameters are also nessesary (CADD damage score threshold, number of permutations, number of samples to be used as controls, etc), but this parameters are asked from user to enter interactivly. 
This script has two modes, interactive (recomanded for few permutations) or you can send it like a job by configuring bash script simmilar as run.tessting.rewas_af.0.005.sh, but with your cluster specifications (recomanded for many permutations).
We would recommand to running it for the first time with one or few permutations to get familliar with nessesary parameters. Then to send it as a job to cluster if nessesary (and you have enough cases compared to controls), and run it for 100 or more perutations. 

Reason to run many permutations is to ensure that you are filtering out frequent population specific SNPs or indels (not present at high frequency at EVS of 1000 genome project) if you have more then enough controls (controls have to be well mached to cases and good representation of population, no much mixture). If you would calculate AF on whole set of controls and filter based on this AF, you would bias your results as you would use same contorls samples for estimating AF and for testing in variant associaition tests.  

For more information how to run this REWAS testing script type at your terminal:
```Shell
Rscript germline_risk_test_v2.R --h
```

After running this script you will have results of the mentioned tests.


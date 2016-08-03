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
As database has restricted access, we recomand to just provide file path and column number in which HUGO gene names are in this file (optional you can print out resuls by setting -v T). This will cause script to get mappings from url link, which is bit slower but is only nessesary to do once.

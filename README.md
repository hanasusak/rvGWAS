# REWAS

---

In this repository are sequence of R scripts to perform REWAS (Rare-variant Exome Wide Association Study). To perform this analysis you will need R and installed nessesary packages (from R console instructions):
```R
install.packages(c("argparse,"ggplot2", "mixOmics","reshape2"))
```

To start you need to perform inital QC check. For that you need at least two input files, multi-sample VCF-like file (with additional annotation, done by eDiVa) and Samples info file (tab separated, where Samples IDs are matching the one in multi-sample file). Additionaly you can add separete mutlis-sample file for InDels and samples IDs correction file (with mappings) in case some sample IDs are not matching in two before mentioned files.

For more information how to run this:
```R
Rscript QC_and_filtering_multisaple_call_v2.R --h
```

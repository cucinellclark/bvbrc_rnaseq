# RNA-Seq Analysis Service

## Overview

The RNA-Seq Analysis Service provides services for aligning, assembling, and testing differential expression on RNA-Seq data. The service provides two recipes for processing RNA-Seq data: 1) Tuxedo, based on the tuxedo suite of tools (i.e., Bowtie, Cufflinks, Cuffdiff); and 2) and HISAT2 for host (human, etc.) reference genomes. The service provides SAM/BAM output for alignment, tab delimited files profiling expression levels, and differential expression test results between conditions. A tutorial for using the RNA-Seq Analysis Service is available here.

The RNA-Seq Service can be accessed from the Services Menu at the top of the BV-BRC website page and via the Command Line Interface (CLI).



## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

This module provides the following application specfication(s):
* [RNASeq](app_specs/RNASeq.md)


## See also

* [RNA-Seq Analysis Service Quick Reference](https://www.bv-brc.org/docs/quick_references/services/rna_seq_analysis_service.html)
* [RNA-Seq Analysis Service](https://www.bv-brc.org/docs/https://bv-brc.org/app/Rnaseq.html)
* [RNA-Seq Analysis Service Tutorial](https://www.bv-brc.org/docs//tutorial/rna_seq/rna_seq.html)



## References

* Kim, D., et al., TopHat2: accurate alignment of transcriptomes in the presence of insertions, deletions and gene fusions. Genome Biol, 2013. 14(4): p. R36.
* Kim, D., B. Langmead, and S.L. Salzberg, HISAT: a fast spliced aligner with low memory requirements. Nat Methods, 2015. 12(4): p. 357-60.
* McClure, R., et al., Computational analysis of bacterial RNA-Seq data. Nucleic Acids Res, 2013. 41(14): p. e140.


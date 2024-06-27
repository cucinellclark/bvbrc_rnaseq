# RNA-Seq Analysis Service

## Overview

The RNA-Seq Analysis Service provides services for aligning, assembling, and testing differential expression on RNA-Seq data. The service provides three recipes for processing RNA-Seq data: 1) Tuxedo, based on the tuxedo suite of tools (i.e., Bowtie, Cufflinks, Cuffdiff) and 2) Bowtie, HTSeq, and DESeq2 for bacterial reference genomes; and 3) and HISAT2, Stringtie, and DESeq2 for host (human, etc.) reference genomes. The service provides SAM/BAM output for alignment, tab delimited files profiling expression levels, and differential expression test results between conditions. A tutorial for using the RNA-Seq Analysis Service is available [here](https://www.bv-brc.org/docs/tutorial/rna_seq/rna_seq.html).

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
* Anders, S., et al., HTSeq—a Python framework to work with high-throughput sequencing data. Bioinformatics, 2014. 31(2): p. 166-169.
* Love, M., et al., Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 2014. 15(550)
* Pertea, M., et al., StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. nature biotechnology, 2015. 33: p. 290-295
* Wang, L., et al., RSeQC: quality control of RNA-seq experiments. Bioinformatics, 2012. 28(16): p. 2184-2185
* Ewels, P., et al., MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 2016. 32(19): p. 3047-3048

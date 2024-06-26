
# Application specification: RNASeq

This is the application specification for service with identifier RNASeq.

The backend script implementing the application is [App-RNASeq.pl](../service-scripts/App-RNASeq.pl).

The raw JSON file for this specification is [RNASeq.json](RNASeq.json).

This service performs the following task:   Align or assemble RNASeq reads into transcripts with normalized expression levels

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| experimental_conditions | Experimental conditions | string  |  |  |
| contrasts | Contrast list | string  |  |  |
| strand_specific | Are the reads in this study strand-specific? | bool  |  | 1 |
| paired_end_libs |  | group  |  |  |
| single_end_libs |  | group  |  |  |
| srr_libs |  | group  |  |  |
| reference_genome_id | Reference genome ID | string  | :heavy_check_mark: |  |
| genome_type | genome type: bacteria or host | enum  | :heavy_check_mark: |  |
| recipe | RNAseq recipe | enum  | :heavy_check_mark: | HTSeq-DESeq |
| host_ftp | Host FTP | string  |  |  |
| output_path | Output Folder | folder  | :heavy_check_mark: |  |
| output_file | File Basename | wsid  | :heavy_check_mark: |  |
| trimming | Trimgalore | bool  |  | 0 |
| unit_test | Unit Test Path | string  |  |  |
| skip_sampling | Skip Sampling | string  |  |  |


#!/usr/bin/python3

#library modules
import sys, os, subprocess, shutil
import concurrent.futures
import glob, json
import tempfile

import pandas as pd
import numpy as np

from bvbrc_api import getSubsystemsDf,getPathwayDf

class DifferentialExpression:

    comparisons = None
    genome = None

    def __init__(self, c):
        print("Creating DifferentialExpression manager")
        self.comparisons = c

    # used to set correct counts matrix delimiter and 
    # to run gene+transcript diffexp if host (just gene for bacteria)
    def set_genome(self, g):
        self.genome = g

    def create_metadata_file(self, sample_list, output_dir):
        meta_file = os.path.join(output_dir,'sample_metadata.tsv')
        output_text_list = ['Sample\tCondition']
        for sample in sample_list:
            text = sample.get_id() + '\t' + sample.get_condition()
            output_text_list.append(text)
        output_text = '\n'.join(output_text_list) + '\n'
        with open(meta_file,'w') as o:
            o.write(output_text)
        return meta_file

    # TODO: don't let colons exist in the condition names on UI
    def run_differential_expression(self,output_prefix): 
        contrast_list = self.comparisons.get_contrast_list() 
        gene_counts = self.genome.get_genome_data(self.genome.get_id()+"_gene_counts")
        genome_type = self.genome.get_genome_type()
        meta_file = self.genome.get_genome_data('sample_metadata_file')
        
        deseq_cmd = ['run_deseq2',gene_counts,meta_file,output_prefix,self.genome.get_genome_data('report_img_path'),genome_type]
        contrast_file_list = []
        for contrast in contrast_list:
            deseq_cmd = deseq_cmd + [contrast]
            cond1 = contrast.split(':')[0]
            cond2 = contrast.split(':')[1]
            diffexp_file = output_prefix + cond1 + '_vs_' + cond2 + '.' + genome_type + '.deseq2.tsv' 
            contrast_file_list.append(diffexp_file)
        try:
            print('Running Command:\n{0}'.format(' '.join(deseq_cmd)))
            subprocess.check_call(deseq_cmd)
            self.genome.add_genome_data('contrast_file_list',contrast_file_list)
        except Exception as e:
            sys.stderr.write('Error running run_deseq2:\n{0}'.format(e))
            return -1
         
        # TODO: run transcript counts
        if self.genome.get_genome_type() == 'host':
            print('implement')
            return
            transcript_counts = self.genome.get_genome_data(self.genome.get_id()+"_transcript_counts")

class GenomeData:

    genome = None

    def __init__(self):
        print("Creating GenomeData manager")
    
    def set_genome(self,g):
        self.genome = g

    def run_queries(self, output_dir, session):
        self.run_subsystems(output_dir, session)
        self.run_pathway(output_dir, session)

    def create_system_figures(self, output_dir):
        superclass_mapping = self.genome.get_genome_data('superclass_mapping')
        pathway_mapping = self.genome.get_genome_data('pathway_mapping')
        genome_counts = self.genome.get_genome_data("tpm")
        metadata = self.genome.get_genome_data('sample_metadata_file')
        superclass_figure = os.path.join(self.genome.get_genome_data('report_img_path'),self.genome.get_id()+"_Superclass_Distribution")
        pathway_figure = os.path.join(self.genome.get_genome_data('report_img_path'),self.genome.get_id()+"_PathwayClass_Distribution")
        superclass_cmd = ["grid_violin_plots",superclass_mapping,genome_counts,metadata,superclass_figure]
        pathway_cmd = ["grid_violin_plots",pathway_mapping,genome_counts,metadata,pathway_figure]

        try:
            print('Running command:\n{0}'.format(' '.join(superclass_cmd)))
            subprocess.check_call(superclass_cmd) 
            self.genome.add_genome_data('superclass_figure',superclass_figure+'.svg')
        except Exception as e:
            sys.stderr.write('Error creating superclass violin plots:\n{0}\n'.format(e))

        try:
            print('Running command:\n{0}'.format(' '.join(pathway_cmd)))
            subprocess.check_call(pathway_cmd) 
            self.genome.add_genome_data('pathway_figure',pathway_figure+'.svg')
        except Exception as e:
            sys.stderr.write('Error creating pathway violin plots:\n{0}\n'.format(e))

    def run_pathway(self, output_dir, session):
        pathway_df = getPathwayDf([self.genome.get_id()], session)
        if not pathway_df is None:
            mapping_table = pathway_df[['patric_id','pathway_class']]
            mapping_output = os.path.join(output_dir,self.genome.get_id()+"_pathway_mapping.tsv")
            mapping_table.to_csv(mapping_output,sep='\t',index=False)
            self.genome.add_genome_data('pathway_mapping',mapping_output)
        else:
            sys.stderr.write('Error, pathway_df is None')
            return -1

    def run_subsystems(self, output_dir, session):
        subsystem_df = getSubsystemsDf([self.genome.get_id()],session)
        if not subsystem_df is None:
            mapping_table = subsystem_df[['patric_id','superclass']]
            mapping_output = os.path.join(output_dir,self.genome.get_id()+"_superclass_mapping.tsv")
            mapping_table.to_csv(mapping_output,sep='\t',index=False)
            self.genome.add_genome_data('superclass_mapping',mapping_output)
        else:
            sys.stderr.write('Error, subsystem_df is None')
            return -1

class Quantify:

    genome = None

    def __init__(self):
        print("Creating Quantify manager") 
    
    def set_genome(self,g):
        self.genome = g

    def run_quantification(self, sample_list, recipe, threads):
        if recipe == 'HTSeq-DESeq':
            htseq_ret = self.run_htseq(sample_list, threads)
            if htseq_ret != 0:
                return htseq_ret
            return self.run_tpmcalc(sample_list, threads)
        elif recipe == 'Host':
            self.run_stringtie(sample_list, threads)
        elif recipe == 'cufflinks':
            self.run_cufflinks(sample_list, threads)            
        else:
            sys.stderr.write("Invalid recipe: {0}".format(recipe))
            return -1 

    def run_tpmcalc(self, sample_list, threads):
        tpm_calc_list = []
        sample_details_list = []
        genome_gtf = self.genome.get_genome_data('gtf')
        for sample in sample_list:
            bam_file = sample.get_sample_data('bam')
            tpm_cmd = ['TPMCalculator','-g',genome_gtf,'-b',bam_file]
            if self.genome.get_genome_type() == 'host':
                tpm_cmd += ['-e']
            if self.genome.get_genome_type() == 'bacteria':
                tpm_cmd += ['-k','gene_name']            
            tpm_output = os.path.join(self.genome.get_sample_path(sample.get_id()),'tpm_counts.tsv') 
            sample_details_list.append([tpm_output,sample])
            tpm_calc_list.append(tpm_cmd) 
        tpm_args_list = list(zip(tpm_calc_list, sample_details_list))
        future_returns = []
        # TPMCalculator outputs to current directory
        if not os.path.exists('TPMCalculator'):
            os.mkdir('TPMCalculator')
        os.chdir('TPMCalculator')
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as pool:
            future_returns = list(pool.map(self.run_tpmcalc_job,tpm_args_list))
        os.chdir('../')
        for f in future_returns:
            if f != 0:
                sys.stderr.write('Error in HTSeq-count: check logs\n')
                sys.stderr.write('{0}\n'.format(future_returns))
                return -1
        return 0 

    def run_tpmcalc_job(self, cmd_details):
        cmd = cmd_details[0]
        sample_details = cmd_details[1]
        sample = sample_details[1]
        sample.add_command("tpmcalc"+"_"+self.genome.get_id(),cmd,"running")
        print('Running command:\n{0}\n'.format(' '.join(cmd)))
        try:
            subprocess.check_call(cmd)
            sample.set_command_status("tpmcalc"+"_"+self.genome.get_id(),"finished")
            output_file = os.path.abspath(sample.get_id()+'_genes.out')
            if not os.path.exists(output_file):
                sys.stderr.write('TPMCalculator output file does not exist:\n{0}\n'.format(output_file))
            sample.add_sample_data(self.genome.get_id()+"_tpm_out",output_file)
        except Exception as e:
            sys.stderr.write('Error running concurrent tpm job:{0}'.format(e))
            sample.set_command_status("tpm"+"_"+self.genome.get_id(),e)
            return -1
        return 0

    def run_htseq(self, sample_list, threads):
        # TODO: add strandedness parameter: -s
        # featurey_type: CDS or Gene
        quant_cmd_list = [] 
        annotation_file = self.genome.get_genome_data('annotation')
        sample_details_list = []
        for sample in sample_list:
            bam_file = sample.get_sample_data('bam')
            quant_cmd = ['htseq-count','-t','gene','-f','bam','-r','pos','-i','ID',bam_file,annotation_file]
            quant_cmd_list.append(quant_cmd)
            sample_dir = self.genome.get_sample_path(sample.get_id())
            sample_output_file = os.path.join(sample_dir,sample.get_id()+'.counts')
            sample_details_list.append([sample_output_file, sample])
        quant_args_list = list(zip(quant_cmd_list,sample_details_list)) 
        future_returns = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as pool:
            future_returns = list(pool.map(self.run_htseq_job,quant_args_list))
        for f in future_returns:
            if f != 0:
                sys.stderr.write('Error in HTSeq-count: check logs\n')
                sys.stderr.write('{0}\n'.format(future_returns))
                return -1
        return 0 
            
    def run_htseq_job(self,cmd_details):
        cmd = cmd_details[0]
        sample_details = cmd_details[1]
        output_file = sample_details[0]
        sample = sample_details[1]
        sample.add_command("htseq"+"_"+self.genome.get_id(),cmd,"running")
        print('Running command:\n{0}\n'.format(' '.join(cmd)))
        try:
            with open(output_file,'w') as o:
                subprocess.check_call(cmd,stdout=o)
            sample.set_command_status("htseq"+"_"+self.genome.get_id(),"finished")
            sample.add_sample_data(self.genome.get_id()+"_gene_counts",output_file)
        except Exception as e:
            sys.stderr.write('Error running concurrent htseq job:{0}'.format(e))
            sample.set_command_status("htseq"+"_"+self.genome.get_id(),e)
            return -1
        return 0
    
    def create_genome_counts_table(self,output_dir,sample_list):
        if self.genome.get_genome_type() == 'bacteria':
            return self.create_genome_counts_table_htseq(output_dir,sample_list)
        elif self.genome.get_genome_type() == 'host':
            return self.create_genome_counts_table_stringtie(output_dir,sample_list)

    # Call prepDE.py script which formats data in the gtf files as a table
    # https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
    def create_genome_counts_table_stringtie(self,output_dir,sample_list):
        # create input file for -i option
        output_text_list = []
        avg_len_list = []
        for sample in sample_list:
            gtf_path = sample.get_sample_data(self.genome.get_id()+'_merged_transcripts')
            text = sample.get_id() + ' ' + gtf_path
            output_text_list.append(text)
            avg_len_list.append(int(sample.get_sample_data('avg_read_length')))
        output_text = '\n'.join(output_text_list) + '\n'
        path_file = os.path.join(output_dir,'sample_transcript_paths.txt')
        with open(path_file,'w') as o:
            o.write(output_text)
        # run prepDE.py 
        gene_matrix_file = os.path.join(output_dir,self.genome.get_id()+"_gene_counts.csv")
        transcript_matrix_file = os.path.join(output_dir,self.genome.get_id()+"_transcript_counts.csv")
        avg_read_length = str(int(np.average(avg_len_list))) 
        # TODO: set path or import or something
        prepde_cmd = ['prepDE','-i',path_file,'-g',gene_matrix_file,'-t',transcript_matrix_file,'l',avg_read_length]
        try:
            print(' '.join(prepde_cmd))
            subprocess.check_call(prepde_cmd)
            self.genome.add_genome_data(self.genome.get_id()+"_gene_counts",gene_matrix_file)
            self.genome.add_genome_data(self.genome.get_id()+"_transcript_counts",transcript_matrix_file)
        except Exception as e:
            sys.stderr.write('Error in prepDE.py: cannot generate genome counts or transcript counts file') 
            sys.exit(-1)
        
    def create_genome_counts_table_htseq(self,output_dir,sample_list):
        genome_df = None
        for sample in sample_list:
            sample_df = pd.read_csv(sample.get_sample_data(self.genome.get_id()+"_gene_counts"),delim_whitespace=True,index_col=0,header=None,names=[sample.get_id()])
            if genome_df is None:
                genome_df = sample_df
            else:
                genome_df = genome_df.join(sample_df)
        output_file = os.path.join(output_dir,self.genome.get_id() + '_counts.tsv')
        genome_df.to_csv(output_file,sep='\t')
        self.genome.add_genome_data(self.genome.get_id()+'_gene_counts',output_file)
        return output_file

    def create_genome_tpm_table(self, output_dir, sample_list):
        if self.genome.get_genome_type() == 'bacteria':
            return self.create_tpm_table_tpmcalculator(output_dir,sample_list)
        if self.genome.get_genome_type() == 'host':
            return self.create_tpm_table_stringtie(output_dir,sample_list)

    def create_tpm_table_tpmcalculator(self, output_dir, sample_list):
        genome_df = None
        for sample in sample_list:
            sample_df = pd.read_csv(sample.get_sample_data(self.genome.get_id()+'_tpm_out'),delim_whitespace=True)
            sample_df = sample_df[['Gene_Id','TPM']]
            sample_df.rename(columns={'Gene_Id':'Gene_Id','TPM':sample.get_id()}, inplace=True)
            if genome_df is None:
                genome_df = sample_df
            else:
                genome_df = genome_df.merge(sample_df,how='outer',on='Gene_Id')
        genome_df = genome_df.fillna(0)
        genome_df.set_index('Gene_Id',inplace=True)
        output_file = os.path.join(output_dir,self.genome.get_id()+"_tpm_counts.tsv")
        genome_df.to_csv(output_file,sep='\t')
        self.genome.add_genome_data('tpm',output_file)

    def create_tpm_table_stringtie(self,output_dir,sample_list):
        genome_df = None
        for sample in sample_list:
            sample_df = pd.read_csv(sample.get_sample_data(self.genome.get_id()+'_merged_gene_counts'),sep='\t')
            sample_df = sample_df.loc[sample_df['Gene ID'] != '.']
            sample_df = sample_df[['Gene ID','TPM']]
            sample_df.rename(columns={'Gene ID':'Gene_Id','TPM':sample.get_id()}, inplace=True)
            if genome_df is None:
                genome_df = sample_df
            else:
                genome_df = genome_df.merge(sample_df,how='outer',on='Gene_Id')
            genome_df = genome_df.fillna(0)
            genome_df.set_index('Gene_Id', inplace=True)
            output_file = os.path.join(output_dir,self.genome.get_id()+"_tpm_counts.tsv")
            genome_df.to_csv(output_file,sep='\t')
            self.genome.add_genome_data('tpm',output_file)

    # TODO: for all sample stuff replace sample.get_id() with self.genome.get_id() 
    def run_stringtie(self, sample_list, threads):
        annotation_file = self.genome.get_genome_data('annotation')
        gtf_list = []
        print('sample_list = [{0}],'.format(sample_list))
        for sample in sample_list:
            bam_file = sample.get_sample_data('bam')
            sample_dir = self.genome.get_sample_path(sample.get_id())

            gene_output = os.path.join(sample_dir,'gene_abund.tab')
            gtf_output = os.path.join(sample_dir,'transcripts.gtf')
            gtf_list.append(gtf_output)
            # not including -e (included below), including restricts novel feature prediction 
            quant_cmd = ['stringtie',bam_file,'-p',str(threads),'-A',gene_output,'-G',annotation_file,'-o',gtf_output]
            if not os.path.exists(gtf_output) or True:
                sample.add_command("stringtie_"+self.genome.get_id(),quant_cmd,"running")
                print('Running command:\n{0}\n'.format(' '.join(quant_cmd))) 
                try:
                    subprocess.check_call(quant_cmd)
                    sample.set_command_status('stringtie_'+self.genome.get_id(),'finished')
                    sample.add_sample_data(self.genome.get_id()+'_gene_counts',gene_output)
                    sample.add_sample_data(self.genome.get_id()+'_transcripts',gtf_output)
                except Exception as e:
                    sys.stderr.write('Error running stringtie:\n{0}\n'.format(e))
                    sample.set_command_status('stringtie_'+self.genome.get_id(),e)
                    return -1
            else:
                sys.stderr.write('{0} already exists: skipping stringtie'.format(gtf_output))
        # merge reconstructed transcriptomes
        merge_file = os.path.join(self.genome.get_genome_dir(),self.genome.get_id()+'.merged.gtf')
        if not os.path.exists(merge_file) or True:
            merge_cmd = ['stringtie','--merge','-G',annotation_file,'-o',merge_file] + gtf_list
            try:
                print('Running command:\n{0}\n'.format(' '.join(merge_cmd)))
                subprocess.check_call(merge_cmd)
            except Exception as e:
                sys.stderr.write('ERROR running stringtie-merge:\n{0}'.format(e))
                return -1
        self.genome.add_genome_data('merged_gtf',merge_file) 

        for sample in sample_list:
            bam_file = sample.get_sample_data('bam')
            sample_dir = self.genome.get_sample_path(sample.get_id())

            gene_output = os.path.join(sample_dir,'merged_gene_abund.tab')
            gtf_output = os.path.join(sample_dir,'merged_transcripts.gtf')
            # -e requires using trancripts found in -G (turns off novel feature prediction)
            quant_cmd = ['stringtie','-e',bam_file,'-p',str(threads),'-A',gene_output,'-G',merge_file,'-o',gtf_output]
            if not os.path.exists(gtf_output) or True:
                sample.add_command('stringtie_merged_'+self.genome.get_id(),quant_cmd,'running')
                print('Running command:\n{0}\n'.format(' '.join(quant_cmd)))
                try:
                    subprocess.check_call(quant_cmd)
                    sample.set_command_status('stringtie_merged_'+self.genome.get_id(),'finished')
                    sample.add_sample_data(self.genome.get_id()+'_merged_transcripts',gtf_output)
                    sample.add_sample_data(self.genome.get_id()+'_merged_gene_counts',gene_output)
                except Exception as e:
                    sys.stderr.write('Error running stringtie-merged:\n{0}\n'.format(e))
                    sample.set_command_status('stringtie_merged_'+self.genome.get_id(),e)
                    return -1
            else:
                sys.stderr.write('{0} already exists: skipping stringtie merged annotation'.format(gtf_output))

    def run_cufflinks(self,sample_list,threads):
        reference = self.genome.get_genome_data('fasta')
        annotation = self.genome.get_genome_data('annotation')
        threads = 8
        for sample in sample_list:
            sample_bam = sample.get_sample_data('bam')
            cufflinks_cmd = ['cufflinks','--quiet','-G',annotation,'-b',reference,'-I','50']
            # Review: Memory mapped location on each system
            # Attempt to copy to /dev/shm. cufflinks seeks a lot in the file.
            # If that fails, try tmp.
            #
            bam_to_use = None

            try:
                tmpfd, bam_tmp = tempfile.mkstemp(prefix="CUFFL.",dir="/dev/shm")
                os.close(tmpfd)
                shutil.copy(sample_bam, bam_tmp)
                print ("Copy succeeded to %s" % (bam_tmp))
                bam_to_use = bam_tmp
            except IOError as err:
                os.unlink(bam_tmp)
                bam_to_use = None
                bam_tmp = None
                sys.stderr.write("Can't copy %s to %s: %s\n" % (sample_bam, bam_tmp, err))

            if bam_to_use == None:
                try:
                    tmpfd, bam_tmp = tempfile.mkstemp(prefix="CUFFL.",dir=".")
                    os.close(tmpfd)
                    shutil.copy(sample_bam, bam_tmp)
                    bam_to_use = bam_tmp


                except IOError as err:
                    os.unlink(bam_tmp)
                    bam_to_use = None
                    bam_tmp = None
                    sys.stderr.write("Can't copy %s to %s: %s\n" % (sample_bam, bam_tmp, err))

            
            if bam_to_use == None:
                sys.stderr.write("Can't copy %s to tmp space\n" % (sample_bam))
                bam_to_use = sample_bam 
                bam_tmp = None

            cufflinks_cmd += [bam_to_use]
            cuff_gtf = sample.get_sample_data('bam').replace('.bam','_transcripts.gtf')
            print("Running command:\n{0}".format(" ".join(cufflinks_cmd)))
            try:
                subprocess.check_call(cufflinks_cmd) 
                sample.add_sample_data('cuff_gtf',cuff_gtf)
            except Exception as e:
                sys.stderr.write('Error running cufflinks on sample {0}:\n{1}\n'.format(sample.get_id(),e))
                return False
            
            if bam_tmp != None:
                sys.stderr.write("remove temp %s\n" % (bam_tmp))
                os.unlink(bam_tmp)
        return True

class Alignment: 

    # TODO: check for genome before running functions
    genome = None
    sampled_reads = 2000
    strand_reads = int(sampled_reads / 2)

    def __init__(self):
        print("Creating Alignment manager")
        
    def set_genome(self,g):
        self.genome = g 

    def run_alignment(self, sample, threads):
        # TODO: if exists, add strand param to alignment
        sample_dir = self.genome.get_sample_path(sample.get_id())
        print("sample_dir = {0}".format(sample_dir))
        reads_list = sample.get_reads_as_list()
        sam_file = os.path.join(sample_dir,sample.get_id()+".sam")
        if self.genome.get_genome_type() == 'bacteria':
            align_cmd = ["bowtie2","-x",self.genome.get_genome_data("bowtie_prefix")] 
        else:
            align_cmd = ["hisat2","-x",self.genome.get_genome_data('hisat_prefix'),'--mp','1,0','--pen-noncansplice','20']
        if sample.sample_type == "paired":
            align_cmd += ["-1",reads_list[0],"-2",reads_list[1]]
        elif sample.sample_type == "single":
            align_cmd += ["-U",reads_list[0]]
        align_cmd += ["-S",sam_file,"-p",str(threads)] 
        sample.add_command("align"+"_"+self.genome.get_id(),align_cmd,"running") 
        align_output_file = sam_file.replace('.sam','.align_stdout')
        print("Running command:\n{0}".format(" ".join(align_cmd))) 
        try:
            # capture stdout for stats later
            with open(align_output_file,'w') as o:
                subprocess.check_call(align_cmd,stderr=o)
            # print captured stdout
            with open(align_output_file,'r') as aof:
                print(aof.read())
            sample.add_sample_data(self.genome.get_id()+'_align_stats',align_output_file)
            sample.set_command_status("align"+"_"+self.genome.get_id(),"finished")
        except Exception as e:
            sys.stderr.write("Sample-alignment encountered an error in Sample {0}:\ncheck error log file".format(sample.get_id()))
            sample.set_command_status("align"+"_"+self.genome.get_id(),e)
            return False
        
        bam_file = self.convert_sam_to_bam(sam_file,threads)
        if bam_file:
            sample.add_sample_data("bam",bam_file) 
            return True
        else:
            sys.stderr.write("Bam file does not exist for Sample {0}:\ncheck error log file".format(sample.get_id()))
            return False

    def run_alignment_stats(self, sample, threads):
        sample_dir = self.genome.get_sample_path(sample.get_id())
        sample_bam = sample.get_sample_data("bam")
        # samtools stats 
        stats_cmd = ["samtools","stats","--threads",str(threads),sample_bam]
        stats_output = os.path.join(sample_dir,sample.get_id()+".samtools_stats")
        sample.add_command("samtools_stats_"+self.genome.get_id(),stats_cmd,"running")
        print("Running command:\n{0}".format(" ".join(stats_cmd)))
        try:
            with open(stats_output,"w") as so:
                subprocess.check_call(stats_cmd,stdout=so)
            sample.set_command_status("samtools_stats_"+self.genome.get_id(),"finished")
        except Exception as e:
            sys.stderr.write("Samtools stats encountered an error in Sample {0}:\ncheck error log file".format(sample.get_id()))
            sample.set_command_status("samtools_stats_"+self.genome.get_id(),e)
        
        sample.add_sample_data('avg_read_length',self.get_average_read_length_per_file(stats_output)) 

        # samstat
        samstat_cmd = ["samstat",sample_bam]
        sample.add_command("samstat_"+self.genome.get_id(),samstat_cmd,"running")
        print("Running command:\n{0}".format(" ".join(samstat_cmd)))
        try:
            subprocess.check_call(samstat_cmd)
            sample.set_command_status("samstat_"+self.genome.get_id(),"finished")
        except Exception as e:
            sys.stderr.write("Samstat encountered an error in Sample {0}:\ncheck error log file".format(sample.get_id()))
            sample.set_command_status("samstat_"+self.genome.get_id(),e)

    #Reads the output from samtools stat and grabs the average read length value
    def get_average_read_length_per_file(self,stats_file):
        with open(stats_file,"r") as sf:
            for line in sf:
                if "average length:" in line: 
                    line = line.strip().split()
                    avg_length = line[-1]
                    return avg_length

    # TODO: incorporate checking the sample alignment results
    def run_sample_alignment(self,sample,threads):
        # TODO: here change sample_dir to the genome directory?
        #   - sample.get_path() to genome.get_sample_path(sample.get_path())
        # sample reads
        sample_dir = sample.get_path()  
        reads = sample.get_reads_as_list()
        sampled_reads_list = []
        readNum = 1
        for r in reads:
            sample_file = r.split(".") # TODO: change where this is output?
            sample_file[len(sample_file)-1] = "sampled"
            sample_file.append("fq")
            sample_file = ".".join(sample_file)
            sampled_reads_list.append(sample_file)
            sample_cmd = ["seqtk","sample","-s","42",r,str(self.sampled_reads)]
            sample.add_command("sample"+str(readNum),sample_cmd,"running")
            print("Running command:\n{0}".format(" ".join(sample_cmd))) 
            try:
                with open(sample_file,"w") as so:
                    subprocess.check_call(sample_cmd,stdout=so)
                sample.set_command_status("sample"+str(readNum),"finished")
            except Exception as e:
                sys.stderr.write("Sampling encountered an error in Sample {0}:\ncheck error log file".format(sample.get_id()))   
                sample.set_command_status("sample"+str(readNum),e)
                return False
            readNum = readNum + 1
        print("sampled reads list:\n{0}".format(sampled_reads_list))
    
        # align sampled reads        
        # TODO: enable for host
        sampled_sam = os.path.join(sample_dir,sample.get_id()+"_sample.sam")
        if self.genome.get_genome_type() == 'bacteria':
            sample_align_cmd = ["bowtie2","-x",self.genome.get_genome_data("bowtie_prefix")]
        else:
            sample_align_cmd = ["hisat2","-x",self.genome.get_genome_data('hisat_prefix'),'--mp','1,0','--pen-noncansplice','20']
        if sample.get_type() == "paired":
            sample_align_cmd += ["-1",sampled_reads_list[0],"-2",sampled_reads_list[1]]
        elif sample.sample_type == "single":
            sample_align_cmd += ["-U",sampled_reads_list[0]]
        sample_align_cmd += ["-S",sampled_sam,"-p",str(threads)] 
        sample.add_command("sample_align",sample_align_cmd,"running") 
        print("Running command:\n{0}".format(" ".join(sample_align_cmd))) 
        try:
            subprocess.check_call(sample_align_cmd)
            sample.set_command_status("sample_align","finished")
        except Exception as e:
            sys.stderr.write("Sample-alignment encountered an error in Sample {0}:\ncheck error log file".format(sample.get_id()))
            sample.set_command_status("sample_align",e)
            return False

        # TODO: condition for assigning strandedness?
        infer_cmd = ["infer_experiment.py","-i",sampled_sam,"-r",self.genome.get_genome_data("bed"),"-s",str(self.strand_reads)]
        infer_file = os.path.join(sample_dir,sample.get_id()+"_strand.infer")
        sample.add_command("infer_strand",infer_cmd,"running")
        print("Running command:\n{0}".format(" ".join(infer_cmd)))
        try:
            with open(infer_file,"w") as o:
                subprocess.check_call(infer_cmd,stdout=o) 
            sample.set_command_status("infer_strand","finished")
            sample.add_sample_data("infer_strand_file",infer_file)
        except Exception as e:
            sys.stderr.write("Infer strand encountered an error in Sample {0}:\ncheck error log file".format(sample.get_id()))
            sample.set_command_status("infer_strand",e)
            return False
        strand = self.infer_strand_from_file(sample.get_sample_data("infer_strand_file"))
        sample.add_sample_data("strand",strand)
        
        # TODO: remove sampled files and sampled sam

    def convert_sam_to_bam(self,sam_file, threads):
        bam_file = sam_file.replace(".sam",".bam")
        print("bam_file = {0}".format(bam_file))
        sam_to_bam_cmd = "samtools view -Su " + sam_file + " | samtools sort -o - - -@ " + str(threads) + " > " + bam_file
        print("Running command:\n{0}".format(sam_to_bam_cmd))
        try:
            subprocess.check_call(sam_to_bam_cmd,shell=True)
        except Exception as e:
            sys.stderr.write("Error in converting sam to bam file:\n{0}\n".format(e))
            return None
        index_cmd = "samtools index " + bam_file
        print("Running command:\n{0}".format(index_cmd))
        try:
            subprocess.check_call(index_cmd,shell=True)
        except Exception as e:
            sys.stderr.write("Error indexing bam file:\n{0}\n".format(e))
            return None
        return bam_file

    def infer_strand_from_file(self,infer_file):
        strand = None
        with open(infer_file,"r") as handle: 
            for x in range(0,3):
                next(handle)
            undetermined = float(next(handle).split(":")[1].strip())
            rf_stranded = float(next(handle).split(":")[1].strip())
            fr_stranded = float(next(handle).split(":")[1].strip())
            if undetermined > fr_stranded and undetermined > rf_stranded:
                strand = "undetermined"
            elif abs(fr_stranded - rf_stranded) <= 0.1: #subjective threshold: should work in most cases
                strand = "undetermined"
            else:
                strand = "FR" if fr_stranded > rf_stranded else "RF"  
        return strand

class DiffExpImport:

    genome = None

    def __init__(self):
        print('Creating differential expression import manager')

    def set_genome(self,g):
        self.genome = g

    def write_gmx_file(self,output_dir):    
        contrast_file_list = self.genome.get_genome_data('contrast_file_list') 
        contrast_list = []
        gene_count_dict = {}
        gene_set = set()
        for contrast_file in contrast_file_list:
            contrast_name = os.path.basename(contrast_file).replace(".tsv","")
            contrast_list.append(contrast_name)
            gene_count_dict[contrast_name] = {}
            with open(contrast_file,"r") as cf:
                next(cf)
                for line in cf:
                    gene,baseMean,log2FC,lfcSE,stat,pvalue,padj = line.strip().split("\t")
                    # strip 'gene-' from identifiers for host
                    gene_set.add(gene.replace("gene-",""))
                    gene_count_dict[contrast_name][gene] = log2FC
        gmx_output = os.path.join(output_dir,'gene_exp.gmx')
        self.genome.add_genome_data('gmx',gmx_output)
        # TODO: rewrite this?
        with open(gmx_output,'w') as o:
            o.write("Gene_ID\t%s\n"%"\t".join(contrast_list))
            for gene in gene_set:
                o.write(gene)
                for contrast in contrast_list:
                    if gene in gene_count_dict[contrast]:
                        o.write("\t%s"%gene_count_dict[contrast][gene])
                    else:
                        o.write("\t0")
                o.write("\n")

    def run_diff_exp_import(self,output_dir,map_args):
        gmx_file = self.genome.get_genome_data('gmx') 
        transform_script = 'expression_transform'
        if os.path.exists(gmx_file):
            experiment_path=os.path.join(output_dir, map_args.d)
            subprocess.call(["mkdir","-p",experiment_path])
            transform_params = {"output_path":experiment_path, "xfile":gmx_file, "xformat":"tsv",\
                    "xsetup":"gene_matrix", "source_id_type":"patric_id",\
                    "data_type":"Transcriptomics", "experiment_title":"RNA-Seq", "experiment_description":"RNA-Seq",\
                    "organism":self.genome.get_id()}
            diffexp_json = self.setup_diffexp_json()
            params_file=os.path.join(output_dir, "diff_exp_params.json")
            # TODO: incorporate check for sstring
            with open(params_file, 'w') as params_handle:
                params_handle.write(json.dumps(transform_params))
            convert_cmd=[transform_script, "--ufile", params_file, "--sstring",map_args.sstring, "--output_path",experiment_path,"--xfile",gmx_file]
            print (" ".join(convert_cmd))
            try:
               subprocess.check_call(convert_cmd)
            except(subprocess.CalledProcessError):
               sys.stderr.write("Running differential expression import failed.\n")
               #subprocess.call(["rm","-rf",experiment_path])
               return
            diffexp_obj_file=os.path.join(output_dir, os.path.basename(map_args.d.lstrip(".")))
            with open(diffexp_obj_file, 'w') as diffexp_job:
                diffexp_job.write(json.dumps(diffexp_json))
            
    def setup_diffexp_json(self):
        #job template for differential expression object
        diffexp_json = json.loads("""                    {
                        "app": {
                            "description": "Parses and transforms users differential expression data",
                            "id": "DifferentialExpression",
                            "label": "Transform expression data",
                            "parameters": [
                                {
                                    "default": null,
                                    "desc": "Comparison values between samples",
                                    "id": "xfile",
                                    "label": "Experiment Data File",
                                    "required": 1,
                                    "type": "wstype",
                                    "wstype": "ExpList"
                                },
                                {
                                    "default": null,
                                    "desc": "Metadata template filled out by the user",
                                    "id": "mfile",
                                    "label": "Metadata File",
                                    "required": 0,
                                    "type": "wstype",
                                    "wstype": "ExpMetadata"
                                },
                                {
                                    "default": null,
                                    "desc": "User information (JSON string)",
                                    "id": "ustring",
                                    "label": "User string",
                                    "required": 1,
                                    "type": "string"
                                },
                                {
                                    "default": null,
                                    "desc": "Path to which the output will be written. Defaults to the directory containing the input data. ",
                                    "id": "output_path",
                                    "label": "Output Folder",
                                    "required": 0,
                                    "type": "folder"
                                },
                                {
                                    "default": null,
                                    "desc": "Basename for the generated output files. Defaults to the basename of the input data.",
                                    "id": "output_file",
                                    "label": "File Basename",
                                    "required": 0,
                                    "type": "wsid"
                                }
                            ],
                            "script": "App-DifferentialExpression"
                        },
                        "elapsed_time": null,
                        "end_time": null,
                        "hostname": "",
                        "id": "",
                        "is_folder": 0,
                        "job_output": "",
                        "output_files": [
                        ],
                        "parameters": {
                            "mfile": "",
                            "output_file": "",
                            "output_path": "",
                            "ustring": "",
                            "xfile": ""
                        },
                        "start_time": ""
                    }
        """)
        return diffexp_json


# Preprocessing: fastqc, trimming, sampled alignment(alignment?), strandedness(alignment?)  
class Preprocess:
    
    def __init__(self):    
        print("Creating Preprocess manager")            

    # TODO: implement multithreaded? 
    def run_fastqc(self, sample): 
        sample_dir = sample.get_path()
        reads = sample.get_reads_as_list()
        fastqc_cmd = ["fastqc","--outdir",sample_dir]
        fastqc_cmd += reads
        sample.add_command("fastqc",fastqc_cmd,"running")  
        print("Running command:\n {0}".format(" ".join(fastqc_cmd))) 
        try:
            subprocess.check_call(fastqc_cmd)
            sample.set_command_status("fastqc","finished")
        except Exception as e:
            sys.stderr.write("FastQC encountered an error in Sample {0}:\ncheck error log file".format(sample.get_id()))   
            sample.set_command_status("fastqc",e)
            return False

    # runs trim_galore on samples
    # TODO: option for rerunning fastqc?
    # TODO: option for adapter specification
    def run_trimming(self, sample, threads):
        sample_dir = sample.get_path()
        reads = sample.get_reads_as_list()
        trimmed_reads = []
        trim_cmd = ["trim_galore","--output_dir",sample_dir,"--cores",str(threads)]
        if sample.get_type() == "paired":
            trimmed_reads.append(os.path.join(sample_dir,os.path.basename(reads[0]).split(".")[0]+"_val_1.fq"))
            trimmed_reads.append(os.path.join(sample_dir,os.path.basename(reads[1]).split(".")[0]+"_val_2.fq"))
            trim_cmd += ["--paired"]
        if sample.get_type() == "single":
            trimmed_reads.append(os.path.join(sample_dir,os.path.basename(reads[0]).split('.')[0]+"_trimmed.fq"))
        trim_cmd += reads
        sample.add_command("trim",trim_cmd,"running")
        print("Running command:\n{0}".format(" ".join(trim_cmd)))
        try:
            subprocess.check_call(trim_cmd)
            sample.set_command_status("trim","finished")
            sample.set_reads_list(trimmed_reads)
            for cutadapt_file in glob.glob('./*cutadapt.log'):
                os.remove(cutadapt_file)  
        except Exception as e:
            sys.stderr.write("Trimming encountered an error in Sample {0}:\ncheck error log file".format(sample.get_id()))
            sample.set_command_status("trim",e)
            return False
        # TODO: remove cutadapt log files: _cutadapt.log
        for r in sample.get_reads_as_list(): 
            if not os.path.exists(r):
                print("{0} does not exist".format(r))


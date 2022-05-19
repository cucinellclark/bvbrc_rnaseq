#!/usr/bin/python3

#library modules
import argparse, sys, os, json
import requests
import subprocess

#pipeline modules
import experiment
import process
# import report
# import process_statistics
from bvbrc_api import authenticateByEnv

#https://www.nature.com/articles/s41598-020-76881-x

# valid recipes
valid_recipes = ['HTSeq-DESeq','cufflinks','Host']

def main(genome_list, experiment_dict, tool_params, output_dir, comparisons, session, map_args):
    # setup folder structure and genome databases
    setup(output_dir, experiment_dict, genome_list)
    diffexp_flag = comparisons.check_diffexp() 
    
    ### process data independently of genomes
    # Fastqc
    preprocess = process.Preprocess()
    for condition in experiment_dict:
        for sample in experiment_dict[condition].get_sample_list():
            preprocess.run_fastqc(sample)

    # Trimming
    # TODO: replace threads with tool_params value
    for condition in experiment_dict:
        for sample in experiment_dict[condition].get_sample_list():
            preprocess.run_trimming(sample, 8)

    ### Sampled align against one genome?
    # TODO: replace threads with tool_params value
    # TODO: assess strandedness with one genome?
    alignment = process.Alignment()
    for genome in genome_list:
        alignment.set_genome(genome)
        for condition in experiment_dict:
            for sample in experiment_dict[condition].get_sample_list():
                alignment.run_sample_alignment(sample, 8)

    ### Align against each genome
    for genome in genome_list:
        alignment.set_genome(genome) 
        for condition in experiment_dict:
            for sample in experiment_dict[condition].get_sample_list():
                # TODO: error checking after alignment?
                alignment.run_alignment(sample, 8)
                alignment.run_alignment_stats(sample, 8)

    # HTSeq(bacteria), Stringtie(host)
    # TODO: some sort of check to make sure everything finished
    # TODO: test host paired
    quantifier = process.Quantify()
    for genome in genome_list:
        quantifier.set_genome(genome)
        quantifier.set_recipe(map_args.recipe)
        sample_list = []
        for condition in experiment_dict:
            samples = experiment_dict[condition].get_sample_list()
            sample_list = sample_list + samples 
        print('sample_list = {0}'.format(sample_list))
        condition_output_list = quantifier.run_quantification(sample_list,8)
        genome_quant_file = quantifier.create_genome_counts_table(output_dir, sample_list)
        genome.add_genome_data('counts_table', genome_quant_file)
        # TODO: test host
        genome_quant_file = quantifier.create_genome_quant_table(output_dir, sample_list)

    # sample_list used in function below
    sample_list = []
    for condition in experiment_dict:
        samples = experiment_dict[condition].get_sample_list() 
        sample_list = sample_list + samples

    # Differential expression 
    if diffexp_flag:
        diff_exp = process.DifferentialExpression(comparisons) 
        meta_file = diff_exp.create_metadata_file(sample_list, output_dir)
        diffexp_import = process.DiffExpImport()
        diffexp_import.set_recipe(map_args.recipe)
        for genome in genome_list:
            genome.add_genome_data('sample_metadata_file',meta_file)
            diff_exp.set_genome(genome)
            output_prefix = os.path.join(output_dir,genome.get_id()+"_")
            diff_exp.run_differential_expression(output_prefix)
            if genome.get_genome_type() == 'bacteria':
                diffexp_import.set_genome(genome)
                diffexp_import.write_gmx_file(output_dir)
                diffexp_import.run_diff_exp_import(output_dir,map_args)

    # Queries: subsystems, kegg
    # output files are used in creating figures
    genome_data = process.GenomeData()
    for genome in genome_list:
        if genome.get_genome_type() == 'bacteria':
            genome_data.set_genome(genome)
            genome_data.run_queries(output_dir,session)
            genome_data.create_system_figures(output_dir)
    
    # TODO: for now assuming one genome for statistics and report
    genome = genome_list[0] 

    # TODO: finish up colleting sample statistics and generate report
    # Get statistics for samples
    # using sample_list created in differential expression section
    # stats = process_statistics.Statistics()
    # stats.add_samples(sample_list) 
    # stats.get_sample_statistics(genome, output_dir, diffexp_flag)

    # report_manager = report.ReportManager()
    # report_manager.create_report(genome, stats.get_stats_dict(), output_dir, diffexp_flag)

    # TODO: Add command output and status 
    # TODO: Add file cleanup
    print("done")

# sets up initial condition, sample, genome folder structure
# folder stucture is: output_dir/condition/sample/genome
# TODO: diffexp object folder
def setup(output_dir, experiment_dict, genome_list):
    # create output directory
    if not os.path.exists(output_dir):
        print("Creating output directory: {0}".format(output_dir))
        os.mkdir(output_dir)
    else:
        print("Output directory alread exists: {0}".format(output_dir))    
    
    # create subfolder list and add paths to relevant class objects 
    subfolder_list = []
    for condition in experiment_dict:
        cond_path = os.path.abspath(os.path.join(output_dir,condition))
        experiment_dict[condition].set_path(cond_path)
        subfolder_list.append(cond_path)
        for sample in experiment_dict[condition].get_sample_list():
            sample_path = os.path.abspath(os.path.join(experiment_dict[condition].get_path(),sample.get_id()))
            sample.set_path(sample_path)
            subfolder_list.append(sample.get_path())
            for genome in genome_list:
                # TODO: if only 1 genome, don't create subfolders?
                genome_path = os.path.join(sample.get_path(),genome.get_id()) 
                genome.create_path_entry(sample.get_id(),genome_path) 
                subfolder_list.append(genome.get_sample_path(sample.get_id()))

    # create subfolders
    print("Creating {0} subfolders:".format(len(subfolder_list)))
    for folder in subfolder_list:
        if not os.path.exists(folder):
            print("Creating subfolder: {0}".format(folder))
            os.mkdir(folder)
        else:
            print("Subfolder already exists: {0}".format(folder))

    # create report images folder
    report_img_folder = os.path.join(output_dir,'report_images/')
    if not os.path.exists(report_img_folder):
        os.mkdir(report_img_folder)
    # TODO: assuming 1 genome
    for genome in genome_list:
        genome.add_genome_data('report_img_path',report_img_folder)

    # TODO: genome_data checking, delete genome if not enough data??? Throw errors???
    # setup genome index objects
    '''
    genome_links_dir = os.path.join(output_dir,"genome_links")
    if not os.path.exists(genome_links_dir):
        print("Creating {0}".format(genome_links_dir))
        os.mkdir(genome_links_dir)
    else:
        print("{0} already exists".format(genome_links_dir))
    for genome in genome_list:
        genome_data_dir = os.path.join(genome_links_dir,genome.get_id())
        if not os.path.exists(genome_data_dir):
            print("Creating {0}".format(genome_data_dir))
            os.mkdir(genome_data_dir)
        else:
            print("{0} already exists".format(genome_data_dir))
        genome.setup_genome_database(genome_data_dir) 
    '''
    for genome in genome_list:
        genome.setup_genome_database()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # TODO: update help
    parser.add_argument('--jfile',
            help='json file for job {"reference_genome_id": "1310806.3", "experimental_conditions":\
                    ["c1_control", "c2_treatment"], "output_file": "rnaseq_baumanii_1505311", \
                    "recipe": "RNA-Rocket", "output_path": "/anwarren@patricbrc.org/home/test",\
                    "paired_end_libs": [{"read1": "/anwarren@patricbrc.org/home/rnaseq_test/MHB_R1.fq.gz",\
                    "read2": "/anwarren@patricbrc.org/home/rnaseq_test/MHB_R2.fq.gz", "condition": 1},\
                    {"read1": "/anwarren@patricbrc.org/home/rnaseq_test/MERO_75_R1.fq.gz",\
                    "read2": "/anwarren@patricbrc.org/home/rnaseq_test/MERO_75_R2.fq.gz", "condition": 2}], "contrasts": [[1, 2]]}', required=True)
    parser.add_argument('-o', help='output directory. defaults to current directory.', required=False, default=None)
    parser.add_argument('-g', help='csv list of directories each containing all genome data: fna, gff, and hisat indices', required=True)
    parser.add_argument('--sstring', help='json server string specifying api {"data_api":"url"}', required=False, default=None)
    parser.add_argument('-p', help='tool parameters', required=False, type=str,default="{}")
    parser.add_argument('-d', help='differential expression folder', required=False, default='.diff_exp')
    # TODO:
    # link to genome files (gff, fa)??? 
    # link to hisat2 index
    # parameters as input json file OR inline
    # unit testing parameters 

    map_args = parser.parse_args()
    print("map_args:\n{0}".format(map_args))

    # set path to use correct version of samtools
    os.environ['PATH'] = "/disks/patric-common/runtime/samtools-1.9/bin:"+os.environ['PATH']
    os.environ['R_LIBS'] = "/disks/patric-common/runtime/lib/R/library:"+os.environ['R_LIBS']
    print('R_LIBS path = {0}'.format(os.environ['R_LIBS']))

    # output directory
    output_dir = map_args.o
    output_dir = os.path.abspath(output_dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # load library
    job_data = None
    try:
        with open(map_args.jfile, 'r') as job_handle:
            job_data = json.load(job_handle)
    except Exception as e:
        print("Error in opening job json file:\n{0}".format(e))
        sys.exit(-1) # Exception: issue in opening job json file
    if not job_data:
        print("job_data is null")
        sys.exit(-1) # Exception: issue with loading job_data

    # Setup session
    s = requests.Session()
    authenticateByEnv(s)

    # load genome ids
    genome_list = []
    #for i in range(0,len(job_data["reference_genome_id"])):
        #genome_list.append(experiment.Genome(job_data["reference_genome_id"][i],job_data["genome_type"][i],s))
    genome_list.append(experiment.Genome(job_data["reference_genome_id"],job_data["genome_type"],s))

    # DOWNLOAD GENOME DATA: remove from perl side

    # Load genome data
    # Note: not making the assumption the genome directories supplied will be in the same order
    #       as the genomes in reference_genome_id
    genome_dir_list = map_args.g.strip().split(',')
    for genome_dir in genome_dir_list:
        if genome_dir.endswith("/"):
            genome_key = os.path.basename(os.path.dirname(genome_dir))
        else:
            genome_key = os.path.basename(genome_dir)
        for genome in genome_list:
            if genome_key == genome.get_id():
                print(os.listdir(genome_dir))
                for f in os.listdir(genome_dir):
                    data_key = None
                    if f.endswith(".fna") or f.endswith(".fa") or f.endswith(".fasta"):
                        data_key = "fasta"
                    elif f.endswith(".gff"): # TODO: other annotation types?
                        data_key = "annotation"
                    elif f.endswith(".ht2.tar"):
                        data_key = "hisat_index"
                    else:
                        continue
                    # TODO: any file linking?
                    genome.set_genome_dir(genome_dir)
                    genome.add_genome_data(data_key,os.path.abspath(os.path.join(genome_dir,f)))

    # sample_list = [] # maybe don't store this, access samples by condition like in original
    experiment_dict = {}
    condition_list = []
    for cond_str in job_data['experimental_conditions']:
        condition_list.append(cond_str)
        new_condition = experiment.Condition(cond_str)
        experiment_dict[cond_str] = new_condition

    #paired_end_libs
    if 'paired_end_libs' in job_data:
        for paired_sample in job_data['paired_end_libs']:
            sample_reads = [paired_sample['read1'],paired_sample['read2']]
            new_sample = experiment.Sample(paired_sample['sample_id'],'paired',sample_reads,None,paired_sample['condition'])
            experiment_dict[new_sample.condition].add_sample(new_sample)

    # single_end_libs
    if 'single_end_libs' in job_data:
        for single_sample in job_data['single_end_libs']:
            sample_read = [single_sample['read']]
            new_sample = experiment.Sample(single_sample['sample_id'],'single',sample_read,None,single_sample['condition'])
            experiment_dict[new_sample.condition].add_sample(new_sample)

    # TODO: test this
    # put sra-fastq files in <output_dir>/SRA_Fastq/
    # put sra-metadata files in <output_dir>/SRA_Meta/
    # set type to paired or single here
    if 'sra_libs' in job_data:
        sra_fastq_dir = os.path.join(output_dir,'SRA_Fastq')
        sra_meta_dir = os.path.join(output_dir,'SRA_Meta')
        if not os.path.exists(sra_fastq_dir):
            os.mkdir(sra_fastq_dir)
        if not os.path.exists(sra_meta_dir):
            os.mkdir(sra_meta_dir)
        for sra_sample in job_data['sra_libs']:
            srr_id = sra_sample['srr_accession'] 
            meta_file = os.path.join(sra_meta_dir,srr_id+'_meta.txt')
            reads_dir = {}
            try:
                subprocess.check_call(['p3-sra','--out',sra_fastq_dir,'--metadata-file',meta_file,'--id',srr_id])
                with open(meta_file,'r') as meta_handle:
                    job_meta = json.load(meta_handle)
                    files = job_meta[0].get('files',[])
                    for i,f in enumerate(files):
                        if f.endswith('_2.fastq'):
                            read2 = os.path.join(sra_fastq_dir,f)
                            reads_dir['read2'] = read2
                        elif f.endswith('_1.fastq'):
                            read1 = os.path.join(sra_fastq_dir,f)
                            reads_dir['read1'] = read1
                        elif f.endswith('.fastq'):
                            read = os.path.join(sra_fastq_dir,f)
                            reads_dir['read'] = read
                missing_read = False
                for read_key in reads_dir:
                    if not os.path.exists(reads_dir[read_key]):
                        sys.stderr.write('Error: fastq file doesn\'t exist:\n{0}\n'.format(reads_dir[read_key]))
                        missing_read = True
                        continue
                # TODO: exit or something?
                if missing_read:
                    continue
                if 'read2' in reads_dir:
                    reads_list = [reads_dir['read1'],reads_dir['read2']]
                    new_sample = experiment.Sample(srr_id,'paired',reads_list,None,sra_sample['condition'])
                    experiment_dict[new_sample.condition].add_sample(new_sample)
                else: #single end
                    reads_list = [reads_dir['read']]
                    new_sample = experiment.Sample(srr_id,'single',reads_list,None,sra_sample['condition'])
                    experiment_dict[new_sample.condition].add_sample(new_sample)
                
            except Exception as e:
                sys.stderr.write('Error in downloading SRA ID {0}:\n{1}\n'.format(srr_id,e)) 
                # TODO: continue or exit??
                sys.exit(-1)

    # TODO: tool parameters
    tool_params = json.loads(map_args.p)
    print('tool_params = {0}'.format(tool_params))
    
    # TODO:
    # finish samples object list
    # load comparisons object: comparison_list

    # TODO:
    #   - check if diffexp is turned on, turn off if need be
    comparisons = experiment.Comparison() 
    for con in job_data["contrasts"]:
         comparisons.add_contrast(con[0],con[1])

    # TODO: Check if job_data contains 'cufflinks' flag: if true, run old pipeline

    # change into output directory
    os.chdir(output_dir)

    # set recipe in map_args
    map_args.recipe = job_data['recipe'] 

    # If not cufflinks, run pipeline
    main(genome_list, experiment_dict, tool_params, output_dir, comparisons, s, map_args)

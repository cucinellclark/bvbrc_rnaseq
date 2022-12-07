#!/usr/bin/python3

# library modules
import subprocess, sys, os

class ReportManager:

    def __init__(self):
        print("Creating ReportManger")

    def run_multiqc(self,output_path):
        debug_multiqc = False
        remove_data_dir = True
        force_overwrite = True
        report_name = os.path.join(output_path,'multiqc_report')
        multiqc_cmd = ["multiqc","--flat","-o",".","-n",report_name,"-t","simple","."]
        if remove_data_dir:
            multiqc_cmd += ["--no-data-dir"]
        if force_overwrite:
            multiqc_cmd += ["-f"]
        if debug_multiqc:
            multiqc_cmd += ["--lint"] 
        try:
            subprocess.check_call(multiqc_cmd)
        except Exception as e:
            sys.stderr.write('Error running multiqc:\n{0}'.format(e)) 
            return -1

    def create_report(self, genome, output_dir, experiment_dict, report_stats, workspace_dir, diffexp_flag):
        report_lines = []
        genome_name = genome.get_genome_name() 
        report_html_header = self.create_html_header(genome_name)
        report_lines.append(report_html_header)
        report_lines.append('<body>')

        report_header = "<header>\n<div class=\"report-info pull-right\">\n<span>Report Date:</span>{0}<br>\n<span>Organism:</span>{1}\n</div>\n</header>".format('DATE',genome.get_genome_name())
        report_lines.append(report_header)

        # Intro: which genome, pipeline, number of samples, conditions, and if differential expression was performed
        report_lines.append('<section>\n<h2>Summary</h2>')
        report_summary = self.create_summary(report_stats,genome)
        report_lines.append(report_summary)
        report_lines.append('</section>')

        # TODO: Link to multiqc report for basic sample summary statistics
        # - how to make link aware of other file while in workspace?
        report_lines.append('<section>\n<h2>MultiQC Report Link</h2>')
        report_lines.append(self.create_multiqc_link(workspace_dir))
        report_lines.append('</section>')

        # Sample IDs and their conditions that have issues
        report_lines.append('<section>\n<h2>Sample Summary Table</h2>')
        report_sample_table = self.create_sample_table(experiment_dict,genome)
        report_lines.append(report_sample_table)
        report_lines.append('</section>')

        # Any figures in report_images
        if genome.get_genome_type() == 'bacteria':
            report_lines.append('<section>\n<h2>Subsystems Distribution</h2>')
            report_subsystems = self.get_subsystem_figure(genome)
            report_lines.append(report_subsystems)
            report_lines.append('</section>')
            report_lines.append('<section>\n<h2>Pathways Distribution</h2>')
            report_pathways = self.get_pathway_figure(genome)
            report_lines.append(report_pathways)
            report_lines.append('</section>')

        # differential expression section if turned on
        # TODO: add when svg API issue is fixed
        #if diffexp_flag:

        # Add references
        report_lines.append('<section>\n<h2>References</h2>')
        reference_lines = self.get_references()
        report_lines.append(reference_lines)
        report_lines.append('</section>')

        report_lines.append('</body>')
        report_lines.append('</html>')
        report_html = '\n'.join(report_lines)
        output_report = os.path.join(output_dir,'bvbrc_rnaseq_report.html')
        with open(output_report,'w') as o:
            o.write(report_html)

    def create_summary(self, report_stats, genome):
        summary_str = f"The BV-BRC RNASeq service was run using the {report_stats['recipe']} pipeline with {report_stats['num_samples']} samples and {report_stats['num_conditions']} conditions." 
        summary_str += f" The reference genome used was {genome.get_genome_name().replace('_',' ')} ({genome.get_id()})."
        return summary_str

    def create_multiqc_link(self, workspace_path):
        #url_base = "https://www.bv-brc.org/workspace"
        url_base = ''
        if workspace_path[-1] != '/':
            workspace_path += '/'
        url = url_base + workspace_path + 'multiqc_report.html'
        link_text = f'<a href=\"multiqc_report.html\" target=\"_parent\">multiqc report link</a>'
        return link_text 

    def get_subsystem_figure(self, genome):
        subsystem_figure_path = genome.get_genome_data('superclass_figure') 
        if subsystem_figure_path and os.path.exists(subsystem_figure_path):
            try:
                with open(subsystem_figure_path,'r') as sfp: 
                    subsystem_figure = sfp.read()
                #subsystem_figure = f'<img src = \"report_images/{genome.get_id()}_Superclass_Distribution.svg\"'
                return subsystem_figure
            except Exception as e:
                sys.stderr.write(f'Error reading file: {subsystem_figure_path}\n')
                return '<p>Error: no subsystem figure found</p>'
        else:
            return '<p>Error: no subsystem figure found</p>'

    def get_pathway_figure(self, genome):
        pathway_figure_path = genome.get_genome_data('pathway_figure') 
        if pathway_figure_path and os.path.exists(pathway_figure_path):
            try:
                with open(pathway_figure_path,'r') as pfp: 
                    pathway_figure = pfp.read()
                #pathway_figure = f'<img src = \"report_images/{genome.get_id()}_Pathway_Distribution.svg\"'
                return pathway_figure
            except Exception as e:
                sys.stderr.write(f'Error reading file: {pathway_figure_path}\n')
                return '<p>Error: no pathway figure found</p>'
        else:
            return '<p>Error: no pathway figure found</p>'
    
    def create_sample_table(self,experiment_dict,genome): 
        table_list = []
        table_list.append("<table class=\"sm-table kv-table center\">")
        table_list.append("<thead class=\"table-header\">")
        table_list.append("<tr>\n<th colspan=\"4\">\n<table-num>Table 1.</table-num>Sample Details\n</th>\n</tr>\n</thead>")
        table_list.append("<tbody>")
        #table_list.append("<tr>\n<td>Condition</td>\n<td>Sample</td>\n<td>Quality</td>\n<td>Alignment</td>\n</tr>")
        table_list.append("<tr>\n<td>Condition</td>\n<td>Sample</td>\n<td>Alignment</td>\n</tr>")
        # TODO: add quality and alignment stats
        # TODO: store stats in sample objects
        for condition in experiment_dict:
            if condition == 'no_condition':
                condition_str = 'None'
            else:
                condition_str = condition
            for sample in experiment_dict[condition].get_sample_list():
                #new_line = f"<tr>\n<td>{condition}</td>\n<td>{sample.get_id()}</td>\n<td>QUALITY</td>\n<td>ALIGNMENT</td>"
                align_file = sample.get_sample_data(genome.get_id()+'_align_stats')
                align_str = 'ALIGNMENT'
                if os.path.exists(align_file):
                    with open(align_file,'r') as af:
                        align_text = af.readlines()
                        align_str = align_text[-1].split(' ')[0]
                new_line = f"<tr>\n<td>{condition_str}</td>\n<td>{sample.get_id()}</td>\n<td>{align_str}</td>"
                table_list.append(new_line)
        table_list.append("</tbody>")
        table_list.append('</table>')
        return '\n'.join(table_list)

    def get_references(self):
        print('generating references')
        # multiqc
        multiqc_ref = "\t\t1.\t MultiQC: Summarize analysis results for multiple tools and samples in a single report Philip Ewels, Måns Magnusson, Sverker Lundin and Max Käller Bioinformatics (2016) doi: 10.1093/bioinformatics/btw35PMID: 27312411 doi: 10.1093/bioinformatics/btw354 PMID: 27312411 "
        # fastqc
        fastqc_ref = "\t\t2.\tAndrews, S. (2010). FastQC:  A Quality Control Tool for High Throughput Sequence Data [Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/"
        # trimgalore
        trimgalore_ref = "\t\t3.\tKrueger, F. (2019, November). Babraham bioinformatics - trim galore! Retrieved December 7, 2022, from https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/ "
        # cutadapt
        cutadapt_ref = "\t\t4.\tMARTIN, Marcel. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, [S.l.], v. 17, n. 1, p. pp. 10-12, may 2011. ISSN 2226-6089. Available at: <https://journal.embnet.org/index.php/embnetjournal/article/view/200>. Date accessed: 07 dec. 2022. doi:https://doi.org/10.14806/ej.17.1.200. "
        # htseq
        htseq_ref = "\t\t5.\tAnders S, Pyl PT, Huber W. HTSeq--a Python framework to work with high-throughput sequencing data. Bioinformatics. 2015 Jan 15;31(2):166-9. doi: 10.1093/bioinformatics/btu638. Epub 2014 Sep 25. PMID: 25260700; PMCID: PMC4287950."
        # stringtie
        stringtie_ref = "\t\t6.\tPertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT, Salzberg SL. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nat Biotechnol. 2015 Mar;33(3):290-5. doi: 10.1038/nbt.3122. Epub 2015 Feb 18. PMID: 25690850; PMCID: PMC4643835."
        # deseq2
        deseq_ref = "\t\t7.\tLove MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8."
        # bowtie/hisat
        bowtie_ref = "\t\t8.\tLangmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. Nat Methods. 2012 Mar 4;9(4):357-9. doi: 10.1038/nmeth.1923. PMID: 22388286; PMCID: PMC3322381."
        hisat_ref = "\t\t9.\tKim, D., Paggi, J.M., Park, C. et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol 37, 907–915 (2019). https://doi.org/10.1038/s41587-019-0201-4"
        # tpmcalculator
        tpmcalc_ref = "\t\t10.\tRoberto Vera Alvarez, Lorinc Sandor Pongor, Leonardo Mariño-Ramírez, David Landsman, TPMCalculator: one-step software to quantify mRNA abundance of genomic features, Bioinformatics, Volume 35, Issue 11, 1 June 2019, Pages 1960–1962, https://doi.org/10.1093/bioinformatics/bty896"
        # cufflinks
        cufflinks_ref = "\t\t11.\tTrapnell, C., Williams, B., Pertea, G. et al. Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and isoform switching during cell differentiation. Nat Biotechnol 28, 511–515 (2010). https://doi.org/10.1038/nbt.1621"
        # seqtk
        seqtk_ref = "\t\t12.\tlh3. (n.d.). Lh3/SEQTK: Toolkit for processing sequences in FASTA/Q Formats. GitHub. Retrieved December 7, 2022, from https://github.com/lh3/seqtk "
        # rseqc
        # bvbrc

    def create_html_header(self,genome_name):
        header = "<!DOCTYPE html><html><head>\n"
        header += "<title>BV-BRC RNASeq Report | {0}</title>\n".format(genome_name)
        header += "<link href=\"https://fonts.googleapis.com/css?family=Work+Sans:300,400,500,600,700,800,900\" rel=\"stylesheet\">\n"
        header += """
            <style>
                body {
                    font-family: 'Work Sans', sans-serif;
                    color: #444;
                }
                header { padding-bottom: .5in; }
                section { margin-bottom: .5in; }

                a {
                    text-decoration: none;
                    color: #0d78ef;
                }
                a:hover { text-decoration: underline; }
                h2 { font-size: 1.45em; }
                h2, h3 { font-weight: 500; }
                h2 small { color: #888; font-size: 0.65em; }
                .pull-left { float: left; }
                .pull-right { float: right; }
                sup { display: inline-block; }

                table-num,
                table-ref,
                fig-num,
                fig-ref {
                    font-weight: 600;
                }

                /* tables */
                table {
                    margin: .1in 0 .5in 0;
                    border-collapse: collapse;
                    page-break-inside: avoid;
                }
                .table-header {
                    color: #fff;
                    background: #196E9C;
                    padding: 5px;
                    text-align: left;
                }
                .table-header th { font-weight: 400; }
                .row-header { border-bottom: 1px solid #196E9C; font-weight: 600; }
                th, td { padding: 5px; }
                th  { border-bottom: 2px solid #196E9C; }
                tr:last-child { border-bottom: 1px solid #196E9C; }

                .kv-table,
                .kv-table-2 {
                    text-align: left;
                }

                .kv-table-2 td:nth-child(2) {
                    border-right: 1px solid rgb(172, 172, 172)
                }

                .kv-table td:first-child,
                .kv-table-2 td:first-child,
                .kv-table-2 td:nth-child(3) {
                    font-weight: 700;
                }

                table td.align-right {
                    text-align: right;
                }

                .kv-table-2 td:nth-child(2) { padding-right: 10px; }
                .kv-table-2 td:nth-child(3) { padding-left: 10px; }
                .lg-table { width: 80%; }
                .med-table { width: 60%; }
                .sm-table { width: 40%; }
                .xs-table {width: 20%; }

                .center { margin-left: auto; margin-right: auto; }

                .logo {
                    width: 2.5in;
                    border-right: 3px solid #777;
                    padding-right: 10px;
                    margin-right: 10px;
                }
                .title {
                    padding: 17px 0 0 0px;
                    font-size: 1.3em;
                    color: #777;
                }
                .report-info {
                    text-align: right;
                    margin-right: 30px;
                    font-size:.8em;
                    color: #777;
                }
                .report-info span {
                    font-weight: 600;
                }

                .main-img {
                    width: 250px;
                    margin-right: .2in;
                }

                .subsystem-fig {
                    width: 8in;
                }

                .circular-view-fig {
                    width: 80%;
                }

                .tree-fig {
                    width: 70%;
                }

                ol.references li {
                    margin-bottom: 1em;
                }

                .clearfix:after {
                    visibility: hidden;
                    display: block;
                    font-size: 0;
                    content: " ";
                    clear: both;
                    height: 0;
                }
                .clearfix { display: inline-block; }

                * html .clearfix { height: 1%; }
                .clearfix { display: block; }
            </style>

            </head>
            """
        return header

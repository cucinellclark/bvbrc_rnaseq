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

    def create_report(self, genome, output_dir, diffexp_flag):
        report_lines = []
        genome_name = genome.get_genome_name() 
        report_html_header = self.create_html_header(genome_name)
        report_lines.append(report_header)
        report_lines.append('<body>')

        # Intro: which genome, pipeline, number of samples, conditions, and if differential expression was performed
        report_header = "<header>\n<div class=\"report-info pull-right\">\n<span>Report Date:</span>{0}<br>\n<span>Organism:</span>{1}\n</div>\n</header>".format('DATE',genome)
        report_lines.append(report_header)

        # TODO: Link to multiqc report for basic sample summary statistics
        # - how to make link aware of other file while in workspace?

        # Sample IDs and their conditions that have issues

        # Any figures in report_images

        report_lines.append('</body>')
        report_lines.append('</html>')
        report_html = '\n'.join(report_lines)
        output_report = os.path.join(output_dir,'bvbrc_rnaseq_report.html')
        with open(output_report,'w') as o:
            o.write(output_report)

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

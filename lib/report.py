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

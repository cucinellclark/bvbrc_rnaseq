#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#file format genematrix *.gmx
import os, sys
from math import log


def main(init_args, output_file):
    output_handle=open(output_file, 'w')
    if(len(init_args)<1):
        print ("Usage cuffdiff_to_genematrix.py  <cuffdiff files>")
        exit(0)
    master_list_genes=set()
    master_list_comparisons=set()
    log_lookup={}
    for input_file in init_args:
        input_handle=open(input_file, 'r')
        lines=input_handle.readlines()
        for line in lines[1:]:
            parts=line.strip().split('\t')
            status='NOT OK'
            if len(parts)==14:
                try:
                    #test_id gene_id gene    locus   sample_1        sample_2        status  value_1 value_2 log2(fold_change)       test_stat       p_value q_value significant
                    (gene_col, sample1, sample2, status, value1, value2, log_change) = (parts[2], parts[4], parts[5], parts[6], float(parts[7]), float(parts[8]), float(parts[9]))
                except ValueError:
                    sys.stderr.write('One of the input files does not match the formatting of a CuffDiff gene differential expression testing file\n')
                    sys.exit()
            else:
                sys.stderr.write('One of the input files does not match the formatting of a CuffDiff gene differential expression testing file\n')
                sys.exit()
            if status != 'OK':
                continue
            gene_ids = []
            if ',' in gene_col:
                gene_ids=gene_col.split(',')
            else:
                gene_ids=[gene_col]
            changed=False
            if value1 == 0:
                value1=0.01
                changed=True
            if value2 == 0:
                value2= 0.01
                changed=True
            if changed:
                log_change= log(value2/value1)/log(2)
            #pandas would be better for this
            for gene_id in gene_ids:
                master_list_genes.add(gene_id)
                comp_id=sample1+' vs '+sample2
                master_list_comparisons.add(comp_id)
                if comp_id not in log_lookup:
                    log_lookup[comp_id]={}
                log_lookup[comp_id][gene_id]=log_change 
        input_handle.close()
    comparisons=list(master_list_comparisons)
    comparisons.sort()
    headers=['Gene ID']+comparisons
    genes=list(master_list_genes)
    genes.sort()
    output_handle.write('\t'.join(headers)+"\n")
    for g in genes:
	    value_list=[]
	    for c in comparisons:
	        try:
		        current_val = str(log_lookup[c][g])
	        except:
		        current_val='NaN'
	        value_list.append(current_val)
	    output_handle.write('\t'.join([g]+value_list)+"\n")
    output_handle.close()

if __name__ == "__main__":
        main(sys.argv[1:])

import csv
import os,sys
import json
from optparse import OptionParser
import glob
import numpy as np

__version__="1.0"
__status__ = "Dev"

###############################
def main():

	in_file1 = "/data/share/project-reza/data/generated/homo_sapiens_gene_expression_bgee.csv"
	in_file2 = "/data/share/project-reza/data/generated/homo_sapiens_gene_expression_TCGA.csv"
	in_file3 = "/data/share/project-reza/data/generated/mus_musculus_gene_expression_bgee.csv"
	
	common_tissues = ["colon","kidney","liver","lung","stomach","urinary bladder"]
	common_genes = ["ENSG00000110955","ENSG00000108518","ENSG00000117984","ENSG00000170315","ENSG00000161016","ENSG00000075624"]
	common_genes += ["ENSG00000198888","ENSG00000197746","ENSG00000167658","ENSG00000166710"]

	common_genes += ["ENSMUSG00000025393","ENSMUSG00000018293","ENSMUSG00000007891","ENSMUSG00000019505","ENSMUSG00000003970","ENSMUSG00000029580"]
	common_genes += ["ENSMUSG00000064341","ENSMUSG00000004207","ENSMUSG00000034994","ENSMUSG00000060802"]	
	
	input_file_list = [in_file1,in_file2,in_file3]
	datasets = ['bgee_human','tcga','bgee_mouse']

	i=-1
	for in_file in input_file_list:
		i+=1
		rowcount = 0
		with open(in_file, 'r') as csvfile:
			csv_reader = csv.reader(csvfile, delimiter=',', quotechar='"')
                	for row in csv_reader:
                        	rowcount +=1
                        	if rowcount ==1:
					tissue_index=[]
                                	for item in common_tissues:
                                        	tissue_index.append(row.index(item))
                                	continue
                        	if row[0] not in common_genes:
                                	continue
                        	else:
                                	fpkm_common_tissues=[]
                                	for item in tissue_index:
                                        	fpkm_common_tissues.append(row[item])
					fpkm_common_tissues = map(float, fpkm_common_tissues)
					if 2*np.mean(fpkm_common_tissues) > max(fpkm_common_tissues):
						flag_mean = 0
					else:
						flag_mean = 1
					if np.std(np.log(fpkm_common_tissues)) < 0.5:
						flag_std = 0
					else:
						flag_std = 1
					print row[0], datasets[i],fpkm_common_tissues,flag_mean,flag_std
 
if __name__ == '__main__':
         main()

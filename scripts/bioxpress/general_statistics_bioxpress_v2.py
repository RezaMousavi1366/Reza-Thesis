import gzip
import os
import glob
import math
import csv

########################################
def main():

 	idmap_file1 = "/data/share/project-reza/data/generated/uberon_tissue.csv"
	in_file = "/data/share/project-reza/data/downloads/bioxpress/BioXpress_interface_overall_final_v2.0.csv"

	rowcount=0
        gene_sign = {}
	uberon_lst = set()
	cancer_lst = set()
        with open(in_file, 'r') as csvfile:
                csvreader = csv.reader(csvfile, delimiter=',', quotechar='|')
                for row in csvreader:
                        rowcount += 1
                        if rowcount ==1:
                                index_significant = row.index("Significant")
                                index_uberon = row.index("UBERON_ID")
				index_cancer = row.index("TCGA Cancer")
				index_gene = row.index("Gene")
                                continue
                        else:
				uberon_id = row[index_uberon].strip()
                        	if uberon_id != "-":
                                	uberon_id = uberon_id.split("&")
                                        for item in uberon_id:
                                        	if item not in uberon_lst:
							print item
                                                	uberon_lst.add(item)
                                if row[0].strip()=='' or row[index_gene].strip()=='':
                                        continue
                                else:
					cancer_type = row[index_cancer].strip()
					gene_id = row[index_gene].strip()
					significant = row[index_significant].strip()
					if gene_id not in gene_sign:
						gene_sign[gene_id]=[]
					gene_sign[gene_id].append(significant)
					if cancer_type != "-":
						if cancer_type not in cancer_lst:
							print cancer_type
							cancer_lst.add(cancer_type)
	count = 0
	for gene_id in gene_sign:
		if all(item=="No" for item in gene_sign[gene_id]):
			count += 1

	print "The number of genes: ", len(gene_sign)
	print "The number of tissues: ", len(uberon_lst)
	print "The number of genes that is not diffrentially expressed is: ", count
	print "The number of cancer types: ", len(cancer_lst)

if __name__ == '__main__':
        main()

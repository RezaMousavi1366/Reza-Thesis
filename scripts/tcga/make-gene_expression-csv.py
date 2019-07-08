import gzip
import os
import glob
import math
import csv

########################################
def mean (expressionValues):
    length = len(expressionValues)
    total_sum = sum(expressionValues)
    average = round(total_sum/length,4)
    return average
########################################
def main():

	idmap_file = "/data/share/project-reza/data/generated/human_protein_id_mapping.csv"
	in_dir = "/data/share/project-reza/data/downloads/TCGA/"
	pattern = in_dir + "/unpacked/*"
	file_lists = glob.glob(pattern)
	out_file = "/data/share/project-reza/data/generated/homo_sapiens_gene_expression_TCGA.csv"

	rowcount=0
 	geneid_canonical = {}
	geneid_status = {}
	geneid_genename = {}
        with open(idmap_file, 'r') as csvfile:
                csvreader = csv.reader(csvfile, delimiter=',', quotechar='|')
                for row in csvreader:
                        rowcount += 1
                        if rowcount ==1:
                                continue
                        else:
                                if row[3]=='':
                                        continue
                                else:
                                        if row[3] not in geneid_canonical:
                                                geneid_canonical[row[3]]=[]
                                                geneid_status[row[3]]=[]
                                                geneid_genename[row[3]]=[]
                                        geneid_canonical[row[3]].append(row[0])
                                        geneid_status[row[3]].append(row[1])
                                        geneid_genename[row[3]].append(row[2])
	
	all_tissues = []
	gene_ids = set()
	data_grid = {}
	for file in file_lists:
    		tissue = file.split("_")[-1].lower()
		if tissue.strip() == "colorectal":
			tissue = "colon"
		if tissue.strip() == "bladder":
			tissue = "urinary bladder"
		if tissue.strip() == "thyroid":
                        tissue = "thyroid gland"
    		if tissue not in all_tissues:
        		all_tissues.append(tissue)
    		file_lists1 = glob.glob(file+"//*")
    		for file1 in file_lists1:
        		if ".txt" in file1:   
            			continue
			else:
            			file_lists2 = glob.glob(file1+"//*")
            			for file2 in file_lists2:
                			with gzip.open(file2,'r') as fin:
                    				for line in fin:
                        				gene_id = line.split('\t')[0].strip().split('.')[0]
							if gene_id in geneid_canonical:
                        					if gene_id not in gene_ids:
                            						gene_ids.add(gene_id)
                        					if gene_id not in data_grid:
                            						data_grid[gene_id] = {}
                        					if tissue not in data_grid[gene_id]:
                            						data_grid[gene_id][tissue] = set()
                        					data_grid[gene_id][tissue].add(line.split('\t')[1].strip())
                        
	for gene_id in gene_ids:
    		for tissue in all_tissues:
        		tissues_per_gene = map(float, data_grid[gene_id][tissue])
        		for index in xrange(0,len(tissues_per_gene)):
            			if tissues_per_gene[index]<1:
                			tissues_per_gene[index]=1
            			tissues_per_gene[index]=math.log(tissues_per_gene[index],2)     
        		data_grid[gene_id][tissue] = mean(tissues_per_gene)

	count = 0
	header_list = ["ensembl_gene_id", "uniprotkb_canonical_ac", "status", "genename"] + all_tissues
	with open(out_file, 'wb') as csvfile:
    		writer = csv.writer(csvfile)
    		writer.writerow(header_list)
    		for gene_id in gene_ids:
        		line=[gene_id]
			line += ["|".join(geneid_canonical[gene_id])]
                        line += ["|".join(geneid_status[gene_id])]
                        line += ["|".join(geneid_genename[gene_id])]
        		for tissue in all_tissues:
            			line.append(data_grid[gene_id][tissue])
                        if (all(i>0.1 for i in line[4:])):
                                writer.writerow(line)
                                count += 1
                        elif (sum(line[4:])>0):
                                count += 1
                        else:
                                continue
	print count

if __name__ == '__main__':
        main()

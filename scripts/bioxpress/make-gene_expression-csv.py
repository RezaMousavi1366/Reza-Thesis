import gzip
import os
import glob
import math
import csv

########################################
def main():

	idmap_file1 = "/data/share/project-reza/data/generated/human_protein_id_mapping.csv"
	idmap_file2 = "/data/share/project-reza/data/generated/uberon_tissue.csv"
	in_file = "/data/share/project-reza/data/downloads/bioxpress/BioXpress_interface_overall_final_v2.0.csv"
	out_file = "/data/share/project-reza/data/generated/homo_sapiens_differential_gene_expresslion_bioxpress.csv"
	
	rowcount=0
        geneid_canonical = {}
        geneid_status = {}
        geneid_genename = {}
        with open(idmap_file1, 'r') as csvfile:
                csvreader = csv.reader(csvfile, delimiter=',', quotechar='|')
                for row in csvreader:
                        rowcount += 1
                        if rowcount ==1:
                                continue
                        else:
                                if row[3]=='':
                                        continue
                                else:
					gene_ids = row[3].split('|')
					for item in gene_ids:
                                        	if item not in geneid_canonical:
                                                	geneid_canonical[item]=[]
                                                	geneid_status[item]=[]
                                                	geneid_genename[item]=[]
                                        	geneid_canonical[item].append(row[0])
                                        	geneid_status[item].append(row[1])
                                        	geneid_genename[item].append(row[2])

	common_tissues = ["colon","kidney","liver","lung","stomach","urinary bladder"]
	uberon_tissue = {}

	uberon_tissue['UBERON:0001155'] = common_tissues[0]
	uberon_tissue['UBERON:0002113'] = common_tissues[1]
	uberon_tissue['UBERON:0002107'] = common_tissues[2]
	uberon_tissue['UBERON:0002048'] = common_tissues[3]
	uberon_tissue['UBERON:0000945'] = common_tissues[4]
	uberon_tissue['UBERON:0001255'] = common_tissues[5]

	rowcount=0
 	data_grid = {}
        with open(in_file, 'r') as csvfile:
        	csvreader = csv.reader(csvfile, delimiter=',', quotechar='|')
                for row in csvreader:
                        rowcount += 1
                        if rowcount ==1:
				index_significant = row.index("Significant")
				index_uberon = row.index("UBERON_ID")
                                continue
                        else:
                                if row[0].strip()=='':
                                        continue
                                else:
                                        uberon_id = row[index_uberon].strip()
					if uberon_id in uberon_tissue:
                                        	if row[0] not in data_grid:
                                        		data_grid[row[0]] = {}
							for item in common_tissues:
								data_grid[row[0]][item]=''
                                        	data_grid[row[0]][uberon_tissue[uberon_id]]=row[index_significant]

	rowcount = 0
        header_list = ["ensembl_gene_id", "uniprotkb_canonical_ac", "status", "genename"] + common_tissues
        with open(out_file, 'wb') as csvfile:
                csvwriter = csv.writer(csvfile)
                csvwriter.writerow(header_list)
                for gene_id in geneid_canonical:
			for item in geneid_canonical[gene_id]:
				if item.split('-')[0] in data_grid:
		                        line = [gene_id]
                		        line += [item]
                        		line += ["|".join(geneid_status[gene_id])]
                        		line += ["|".join(geneid_genename[gene_id])]
                        		for tissue in common_tissues:
                                		line.append(data_grid[item.split('-')[0]][tissue])
					csvwriter.writerow(line)

if __name__ == '__main__':
        main()

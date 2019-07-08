import csv
import collections
import math

__version__="1.0"
__status__ = "Dev"

############################################
def mean(expression_values):
	length = len(expression_values)
	if length == 0:
		average=float('nan')
	else:
		total_sum = sum(expression_values)
		average = round(total_sum/length,4)
	return average
############################################
def main():
	
	idmap_file = "/data/share/project-reza/data/generated/human_protein_id_mapping.csv"
	in_file1 = "/data/share/project-reza/data/downloads/bgee_homo_sapiens/homo_sapiens_RNA-Seq_FPKM_GSE30611_Bodymap.tsv"
	in_file2 = "/data/share/project-reza/data/downloads/bgee_homo_sapiens/homo_sapiens_RNA-Seq_FPKM_GSE30352_Brawand.tsv"
	in_file3 = "/data/share/project-reza/data/downloads/bgee_homo_sapiens/homo_sapiens_gene_expression_Fagerberg.csv"
	out_file = "/data/share/project-reza/data/generated/homo_sapiens_gene_expression_bgee.csv"

	tissues_names_Bodymap = ["adipose tissue", "adrenal gland", "kidney", "brain", "colon", "ovary", "heart", "leukocyte", "liver", "lung"]
	tissues_names_Bodymap += ["lymph node", "prostate", "skeletal muscle tissue", "testis", "mammary gland", "thyroid gland"]
	tissues_names_Brawand = ["brain", "heart","kidney","liver", "testis"]
	tissues_names_Fagerberg = ["colon","kidney","liver","pancreas","lung","prostate","brain","stomach","spleen","lymph node","appendix"]
	tissues_names_Fagerberg += ["small intestine","adrenal gland","duodenum","adipose tissue","endometrium","placenta","testis","gallbladder"]
	tissues_names_Fagerberg += ["urinary bladder","thyroid gland","esophagus","heart","skin","ovary","bone marrow","salivary gland"]
	all_tissues = list(set(tissues_names_Bodymap+tissues_names_Brawand+tissues_names_Fagerberg))

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

        rowcount = 0
        gene_ids = set()
	data_grid = {}
        with open(in_file1) as tsvfile:
    		tsvreader = csv.reader(tsvfile, delimiter="\t", quotechar='"')
    		for row in tsvreader:
        		rowcount += 1
        		if rowcount == 1:
            			index_gene_id = row.index("Gene ID")
            			index_tissue = row.index("Anatomical entity name")
            			index_FPKM = row.index("FPKM")
        		else:
            			gene_id = row[index_gene_id]
				if gene_id in geneid_canonical:
					if gene_id not in data_grid:
						gene_ids.add(gene_id)
						data_grid[gene_id] = {}
						for tissue in all_tissues:
							data_grid[gene_id][tissue] = []
					tissue_name = row[index_tissue]
					if tissue_name.strip() == "multi-cellular organism":
						continue
					elif tissue_name.strip() == "adult mammalian kidney":
						data_grid[gene_id]["kidney"].append(row[index_FPKM])
					elif tissue_name.strip() == "thoracic mammary gland":
						data_grid[gene_id]["mammary gland"].append(row[index_FPKM])
					elif tissue_name.strip() == "prostate gland":
						data_grid[gene_id]["prostate"].append(row[index_FPKM])
					elif tissue_name.strip() == "female gonad":
						data_grid[gene_id]["ovary"].append(row[index_FPKM])
					else:
						data_grid[gene_id][tissue_name].append(row[index_FPKM])

	rowcount = 0
	with open(in_file2) as tsvfile:
                tsvreader = csv.reader(tsvfile, delimiter="\t", quotechar='"')
                for row in tsvreader:
                        rowcount += 1
                        if rowcount == 1:
                                index_gene_id = row.index("Gene ID")
                                index_tissue = row.index("Anatomical entity name")
                                index_FPKM = row.index("FPKM")
                        else:
                                gene_id = row[index_gene_id]
				if gene_id in geneid_canonical:
                                	if gene_id not in data_grid:
                                        	gene_ids.add(gene_id)
                                        	data_grid[gene_id] = {}
                                        	for tissue in all_tissues:
                                                	data_grid[gene_id][tissue] = []
                                	tissue_name = row[index_tissue]
					if tissue_name.strip() == "frontal cortex":
						data_grid[gene_id]["brain"].append(row[index_FPKM])
					elif tissue_name.strip() == "prefrontal cortex":
						data_grid[gene_id]["brain"].append(row[index_FPKM])
					elif tissue_name.strip() == "temporal lobe":
						data_grid[gene_id]["brain"].append(row[index_FPKM])
					elif tissue_name.strip() == 'cerebellum':
						data_grid[gene_id]["brain"].append(row[index_FPKM])
					elif tissue_name.strip() == "adult mammalian kidney":
                                                data_grid[gene_id]["kidney"].append(row[index_FPKM])
					else:
						data_grid[gene_id][tissue_name].append(row[index_FPKM])
	
        for gene_id in gene_ids:
                for tissue in all_tissues:
                        tissues_per_genes = map(float, data_grid[gene_id][tissue])
			for index in xrange(0,len(tissues_per_genes)):
				if tissues_per_genes[index]<1:
                                	tissues_per_genes[index]=1
                                tissues_per_genes[index]=math.log(tissues_per_genes[index],2)
                        data_grid[gene_id][tissue] = mean(tissues_per_genes)
	
	with open(in_file3, 'r') as csvfile:
		csvreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        	for row in csvreader:
            		rowcount +=1
            		if rowcount ==1:
                		continue
            		else:
				gene_id = row[0]
                		if gene_id in geneid_canonical:
					if gene_id not in data_grid:
                                		gene_ids.add(gene_id)
                                                data_grid[gene_id] = {}
                                                for tissue in all_tissues:
                                                        data_grid[gene_id][tissue] = []
                    			for idx, val in enumerate(row[1:]):
						tissue_per_gene = data_grid[gene_id][tissues_names_Fagerberg[idx]]
                                                if isinstance(tissue_per_gene, float) and math.isnan(tissue_per_gene)==False:
                                                        data_grid[gene_id][tissues_names_Fagerberg[idx]] = mean([tissue_per_gene,float(val)])
                                                else:
                                                        data_grid[gene_id][tissues_names_Fagerberg[idx]] = round(float(val),4)

	count = 0
	header_list = ["ensembl_gene_id", "uniprotkb_canonical_ac", "status", "genename"] + all_tissues
	with open(out_file, 'wb') as csvfile:
    		csvwriter = csv.writer(csvfile)
		csvwriter.writerow(header_list)
		for gene_id in gene_ids:
			line = [gene_id]
                        line += ["|".join(geneid_canonical[gene_id])]
                        line += ["|".join(geneid_status[gene_id])]
                        line += ["|".join(geneid_genename[gene_id])]
			for tissue in all_tissues:
				line.append(data_grid[gene_id][tissue])
			if (all(i>0.1 for i in line[4:])):
            			csvwriter.writerow(line)
				count += 1
			elif (sum(line[4:])>0):
				count += 1
			else:
				continue
	print count

if __name__ == '__main__':
        main()

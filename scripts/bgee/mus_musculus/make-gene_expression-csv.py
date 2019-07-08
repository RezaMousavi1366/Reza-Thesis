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
	
	idmap_file = "/data/share/project-reza/data/generated/mouse_protein_id_mapping.csv"
	in_file1 = "/data/share/project-reza/data/downloads/bgee_mus_musculus/mus_musculus_RNA-Seq_FPKM_GSE41637_merkin.tsv"
	in_file2 = "/data/share/project-reza/data/downloads/bgee_mus_musculus/mus_musculus_RNA-Seq_FPKM_GSE30352_brawand.tsv"
	in_file3 = "/data/share/project-reza/data/downloads/bgee_mus_musculus/mus_musculus_RNA-Seq_FPKM_GSE30617_keane.tsv"
	in_file4 = "/data/share/project-reza/data/downloads/bgee_mus_musculus/mus_musculus_gene_expression_encode.csv"
	out_file = "/data/share/project-reza/data/generated/mus_musculus_gene_expression_bgee.csv"

	tissues_names_merkin = ["brain","heart","kidney","testis","liver","lung","skeletal muscle tissue","spleen","colon"]
	tissues_names_brawand = ["brain", "heart","kidney","liver", "testis"]
	tissues_names_keane = ["heart","brain","liver","lung","spleen","thymus"]
	tissues_names_encode = ["brain","heart","kidney","liver","lung","placenta","small intestine","spleen","testis","thymus","adrenal gland"]
	tissues_names_encode += ["urinary bladder","colon","duodenum","brain","genitalfatpad","large intestine","mammary gland","ovary"]
	tissues_names_encode += ["subcfatpad","stomach"]
	all_tissues = list(set(tissues_names_merkin+tissues_names_brawand+tissues_names_keane+tissues_names_encode))
	
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
                                        gene_ids = row[3].split('|')
                                        for item in gene_ids:
                				if item not in geneid_canonical:
                    					geneid_canonical[item]=[]
                    					geneid_status[item]=[]
                    					geneid_genename[item]=[]
                				geneid_canonical[item].append(row[0])
                				geneid_status[item].append(row[1])
                				geneid_genename[item].append(row[2])

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
					if tissue_name.strip() == "adult mammalian kidney":
						data_grid[gene_id]["kidney"].append(row[index_FPKM])
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
					if tissue_name.strip() == 'cerebellum':
                                                data_grid[gene_id]["brain"].append(row[index_FPKM])
					elif tissue_name.strip() == "adult mammalian kidney":
                                                data_grid[gene_id]["kidney"].append(row[index_FPKM])
					else:
						data_grid[gene_id][tissue_name].append(row[index_FPKM])
	
	rowcount = 0
        with open(in_file3) as tsvfile:
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
                                        if tissue_name.strip() == "Ammon's horn":
                                                data_grid[gene_id]["brain"].append(row[index_FPKM])
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

	with open(in_file4, 'r') as csvfile:
		csvreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        	for row in csvreader:
            		rowcount += 1
            		if rowcount == 1:
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
						tissue_per_gene = data_grid[gene_id][tissues_names_encode[idx]]
						if isinstance(tissue_per_gene, float) and math.isnan(tissue_per_gene)==False:
							data_grid[gene_id][tissues_names_encode[idx]] = mean([tissue_per_gene,float(val)])
						else:
							data_grid[gene_id][tissues_names_encode[idx]] = round(float(val),4)

	'''for gene_id in gene_ids:
    		for tissue in all_tissues:
        		tissues_per_genes = data_grid[gene_id][tissue]
			if len(tissue_per_genes)>1:
				tissues_per_genes = map(float,tissues_per_genes)
        			data_grid[gene_id][tissue] = mean(tissues_per_genes)'''
	
	count=0
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

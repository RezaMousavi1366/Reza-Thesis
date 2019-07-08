import csv
import os,sys
import json
from optparse import OptionParser
import glob

__version__="1.0"
__status__ = "Dev"

###############################
def load_mapping(map_dict,in_file,dataset,common_tissues,mouse_human_homologous):
	map_dict[dataset] = {}
	tissue_index = []
	rowcount = 0
	with open(in_file, 'r') as csvfile:
        	csv_reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        	for row in csv_reader:
                	rowcount +=1
                	if rowcount ==1:
    				for item in common_tissues:
					tissue_index.append(row.index(item))
				continue
			if row[1]=="":
				continue
			else:
				fpkm_common_tissues=[]
				for item in tissue_index:
					fpkm_common_tissues.append(row[item])
				type_gex_level = specificity_classification(fpkm_common_tissues)
				if type_gex_level == "not expressed":
					continue
				else:
					if dataset == "bgee_mouse":
						if row[0] in mouse_human_homologous:
							human_gene_ids = mouse_human_homologous[row[0]]
							for gene_id in human_gene_ids:
								if gene_id not in map_dict[dataset]:
									map_dict[dataset][gene_id]={}
								if row[0] not in map_dict[dataset][gene_id]:
									map_dict[dataset][gene_id][row[0]]=""
								map_dict[dataset][gene_id][row[0]] = row[1:4]+[type_gex_level]
					else:			
						if row[0] not in map_dict[dataset]:
                                        		map_dict[dataset][row[0]]={}
							map_dict[dataset][row[0]][""] = ""
						map_dict[dataset][row[0]][""] = row[1:4]+[type_gex_level]
###############################
def specificity_classification(fpkm_common_tissues):
	len_tissue = len(fpkm_common_tissues)
	fpkm_common_tissues = map(float, fpkm_common_tissues)
	if all(i == 0 for i in fpkm_common_tissues):
        	type_gex_level = "not expressed"
        elif all(i < 1 for i in fpkm_common_tissues):
                type_gex_level = "not_detected"
	else:
		type_gex_level = []
		max_fpkm_values = max(fpkm_common_tissues)
		max_index_fpkm_values = fpkm_common_tissues.index(max_fpkm_values)
		temp_list = fpkm_common_tissues[0:max_index_fpkm_values]+fpkm_common_tissues[max_index_fpkm_values+1:]
		if all(i*50 < max_fpkm_values for i in temp_list):
			type_gex_level.append("tissue_specific")
		if all(i*5 < max_fpkm_values for i in temp_list):
			type_gex_level.append("tissue_enriched")
		detected_tissues = [item for item in fpkm_common_tissues if item >= 1]
		if len(detected_tissues) < len_tissue:
			if all(i > 6 for i in detected_tissues):
                                type_gex_level.append("mixed_high")
                        else:
                                type_gex_level.append("mixed_low")	
		else:
			if all(i > 6 for i in detected_tissues):
				type_gex_level.append("expressed_in_all_high")
			else:
				type_gex_level.append("expressed_in_all_low")
		type_gex_level = "|".join(type_gex_level)
	return type_gex_level
###############################
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2] 
    return lst3
###############################
def main():
	idmap_file = "/data/share/project-reza/data/generated/human_to_mouse_homologs.csv"
	in_file1 = "/data/share/project-reza/data/generated/homo_sapiens_gene_expression_bgee.csv"
	in_file2 = "/data/share/project-reza/data/generated/homo_sapiens_gene_expression_TCGA.csv"
	in_file3 = "/data/share/project-reza/data/generated/mus_musculus_gene_expression_bgee.csv"
	out_file = "/data/share/project-reza/data/generated/bgee_tcga_gene_expression.csv"

	#common_tissues = ["colon","kidney","liver","lung","prostate","stomach","thyroid gland","esophagus","urinary bladder"]
	common_tissues = ["colon","kidney","liver","lung","stomach","urinary bladder"]
	datasets = ['bgee_human','tcga','bgee_mouse']
	input_file_list = [in_file1,in_file2,in_file3]

	rowcount=0
	mouse_human_homologous = {}
        with open(idmap_file, 'r') as csvfile:
                csvreader = csv.reader(csvfile, delimiter=',', quotechar='|')
                for row in csvreader:
                        rowcount += 1
                        if rowcount ==1:
                                continue
                        else:
                                if row[1]=='':
                                        continue
                                else:
                                        if row[1].replace('"','') not in mouse_human_homologous:
                                        	mouse_human_homologous[row[1].replace('"','')]=[]
                                        mouse_human_homologous[row[1].replace('"','')].append(row[0].replace('"',''))
   
	map_dict={}
	for idx, val in enumerate(datasets):
		load_mapping(map_dict,input_file_list[idx],val,common_tissues,mouse_human_homologous)

	count = 0
	temp = []
	for id1 in map_dict['bgee_mouse']:
		for id2 in map_dict['bgee_mouse'][id1]:
 			if "in_all_high" in map_dict['bgee_mouse'][id1][id2][-1]:
				if id1 not in temp:
					temp.append(id1)
					count += 1
				else:
					print id1
	print count
	print len(id1)
	data_grid={}
	gene_expression_type = []
	for dataset in datasets:
		data_grid[dataset]={}
		for gene_id in map_dict[dataset]:
			for gene_type in map_dict[dataset][gene_id].values():
				list_gene_expression_type = gene_type[-1].split("|")
				for item in list_gene_expression_type:
					if  item not in gene_expression_type:
						gene_expression_type.append(item)
					if item not in data_grid[dataset]:
						data_grid[dataset][item]=[]
					data_grid[dataset][item].append(gene_id)

	print_gene_list = set()
	header_list = ["human_ensembl_gene_id","dataset","mouse_ensembl_gene_id","uniprotkb_canonical_ac","status","gene_name","category"]
	with open(out_file, 'wb') as csvfile:
        	writer = csv.writer(csvfile)
        	writer.writerow(header_list)
		for element in gene_expression_type:
			rowcount = 0
			common_genes = []
			for dataset in datasets:
				rowcount += 1
				if rowcount == 1:
					if element not in data_grid[dataset]:
						lst1_genes = []
					else:
						lst1_genes = data_grid[dataset][element]
					print dataset+"\t"+element+"\t"+str(len(lst1_genes))
				elif rowcount == 2:
					if element not in data_grid[dataset]:
						lst2_genes = []	
					else:
						lst2_genes = data_grid[dataset][element]
					print dataset+"\t"+element+"\t"+str(len(lst2_genes))
					common_genes = [value for value in lst1_genes if value in lst2_genes]
					lst1_genes = []
				else:
					if element not in data_grid[dataset]:
						lst2_genes = []
					else:
						lst2_genes = data_grid[dataset][element]
					print dataset+"\t"+element+"\t"+str(len(lst2_genes))
					common_genes = [value for value in lst2_genes if value in common_genes]
				if len(lst1_genes) != 0:
					gene_ids = lst1_genes
				else:
					gene_ids = lst2_genes
					lst2_genes = []
				for gene_id in gene_ids:
					if dataset == "bgee_mouse":
						if gene_id in print_gene_list:
							continue
						else:
							print_gene_list.add(gene_id)
					for item in map_dict[dataset][gene_id].keys():
						row = [gene_id]
						row += [dataset]
						row += [item]
						row += map_dict[dataset][gene_id][item]
						writer.writerow(row)		
			len_common_genes = len(common_genes)
			print '-'.join(datasets)+"\t"+element+"\t"+str(len_common_genes)
	
if __name__ == '__main__':
         main()

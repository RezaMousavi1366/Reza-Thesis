import csv

import collections
import math

__version__="1.0"
__status__ = "Dev"

###############################
def average(lst): 
    return sum(lst) / len(lst)

###############################
def specificity_classification(fpkm_tissues):

        len_tissue = len(fpkm_tissues)
        fpkm_tissues = map(float, fpkm_tissues)
	fpkm_tissues.sort()
        if len([i for i in fpkm_tissues if i==0])>0:
                type_gex_level = "not expressed in all tissues"
        elif all(i < 1 for i in fpkm_tissues):
                type_gex_level = "not_detected"
        else:
                type_gex_level = []
                max_fpkm_values = max(fpkm_tissues)
                max_index_fpkm_values = fpkm_tissues.index(max_fpkm_values)
                temp_list = fpkm_tissues[0:max_index_fpkm_values]+fpkm_tissues[max_index_fpkm_values+1:]
                if all(i*50 < max_fpkm_values for i in temp_list):
                        type_gex_level.append("tissue_specific")
                if all(i*5 < max_fpkm_values for i in temp_list):
                        type_gex_level.append("tissue_enriched")
		for j in range(2,8):
			small_group = average(fpkm_tissues[0:j])
			large_group = average(fpkm_tissues[j:])
			if small_group > (large_group*5):
				type_gex_level.append("group_enriched")
				break			
                detected_tissues = [item for item in fpkm_tissues if item >= 1]
                if len(detected_tissues) < len_tissue:
                        if all(i > 10 for i in detected_tissues):
                                type_gex_level.append("mixed_high")
                        else:
                                type_gex_level.append("mixed_low")
                else:
                        if all(i > 10 for i in detected_tissues):
                                type_gex_level.append("expressed_in_all_high")
                        else:
                                type_gex_level.append("expressed_in_all_low")
                type_gex_level = "|".join(type_gex_level)
        return type_gex_level
###################################
def main():

	in_file = "/data/share/project-reza/data/generated/homo_sapiens_gene_expression_bgee.csv"
	out_file = "/data/share/project-reza/data/generated/homo_sapiens_gene_expression_classification_bgee.csv"
	
	data_grid = {}
	count = {}
	rowcount = 0
	with open(in_file, 'r') as csvfile:
		csvreader = csv.reader(csvfile, delimiter=',')
        	for row in csvreader:
                	rowcount +=1
                	if rowcount ==1:
				header_list = row + ["category"] 
                                continue
                        else:
				fpkm_values_per_gene = row[4:]
				type_gex_level = specificity_classification(fpkm_values_per_gene)
				gex_level_list = type_gex_level.split("|")
				for item in gex_level_list:
					if item not in count:
						count[item] = 0
					count[item] += 1
				if type_gex_level == "not expressed in all tissues":
                                	continue
                                else:
					if row[0] not in data_grid:
                                		data_grid[row[0]]=row[1:]
                        		data_grid[row[0]].append(type_gex_level)

	for gex_type in count:
		print gex_type+'\t'+str(count[gex_type])
					
        with open(out_file, 'wb') as csvfile:
                csvwriter = csv.writer(csvfile)
                csvwriter.writerow(header_list)
                for gene_id in data_grid:
                        line = [gene_id]
			line += data_grid[gene_id]
                        csvwriter.writerow(line)

if __name__ == '__main__':
        main()

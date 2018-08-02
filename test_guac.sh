#!/usr/bin/bash

data_file='/storage/cylin/home/cl6/projects/a9_rna/data_tables/HMEC_A9_6h_RNA_FASTQ_TABLE.txt'
output_folder='/storage/cylin/home/cl6/projects/a9_rna/guac_out'

gene_dict_path='/storage/cylin/home/cl6/projects/a9_rna/guac_out/hg38_HMEC_A9_gene_gtf.pkl'

#for non overlapping sum159 genes
#gene_list_path='/storage/cylin/home/cl6/projects/a9_rna/tables/HG38_SUM159_FILTERED_GENES_UNIQUE.txt'

#for test genes
gene_list_path='/storage/cylin/home/cl6/projects/a9_rna/tables/HG38_TEST_GENES.txt'

#with a gtf
#python ./1_calculate_pr.py -d $data_file -g hg38 -o $output_folder -title HMEC_A9 -gtf $gene_dict_path -genes $gene_list_path -names A9_6h_Tam_1,A9_6h_Unt_1

#without a gtf
python ./1_calculate_pr.py -d $data_file -g hg38 -o $output_folder -title HMEC_A9 -genes $gene_list_path -names A9_6h_Tam_1,A9_6h_Unt_1

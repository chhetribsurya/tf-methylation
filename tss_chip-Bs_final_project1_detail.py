import numpy as np
import pandas as pd
import pybedtools
import re
import os

#Test_file:with open("/home/surya/head.txt", "r") as file:

input_file_1 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/Pol2_rep1_peak_chr1.bed"
input_file_2 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/Pol2_rep2_peak_chr1.bed"
input_file_3 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/MACS_SL60583_SL60585_10_peaks.bed_chr1.cov"
input_file_4 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/MACS_SL63653_SL63655_10_peaks.bed_chr1.cov"


def intersect_all(*beds, **kwargs):

	#convert args tuple to a list, so as to have the flexibility or
	#mutability advantage of modifying the list later on, if needed:
	bedlist = list(beds)
	x = pybedtools.BedTool(bedlist[0])

	for bed in bedlist[1:]:
		x = x.intersect(bed, **kwargs)
	return x

cleaned_peak_file_1 = intersect_all(input_file_1, input_file_2, input_file_3, input_file_4, wa = True, u = True)
#cleaned_peak_file_1.saveas('/home/surya/surya/for_gata_parse/pybed_overlap.bed')


# You could use functools.reduce to iteratively apply pd.merge to each of the DataFrames:
# result = functools.reduce(merge, dfs)

# This is equivalent to:

# result = dfs[0]
# for df in dfs[1:]:
#     result = merge(result, df)

# - See more at: http://www.mzan.com/article/34338831-pandas-merge-multiple-dataframes-and-control-column-names.shtml#sthash.nx3sNvbM.dpuf


###########

cpg_file_1 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/k562_wgbs_chr1.cov"
cpg_file_2 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/SL60583_GATA2_rep1.cov_chr1.cov"
cpg_file_3 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/SL63653_GATA2_rep2.cov_chr1.cov"
#cpg_file_4 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/SL60583_GATA2_rep1.cov_chr1.cov"
header=["chrom","start","end", "perc_meth", "meth", "unmeth"]

df_1 = pd.read_csv(cpg_file_1, sep="\t", names = header)
df_2 = pd.read_csv(cpg_file_2, sep="\t", names = header)
df_3 = pd.read_csv(cpg_file_3, sep="\t", names = header)
#df_4 = pd.read_csv(cpg_file_4, sep="\t", header=None, names=["chrom","start","end", "meth", "unmeth"])

dfs = [df_1, df_2, df_3]
df_final = reduce(lambda left_df,right_df: pd.merge(left_df,right_df, how='inner', on=['chrom', 'start', 'end']), dfs)

############


#or use the following function for the merging of the cpgs file:
def merge_all(*cpgs, **kwargs):

	#convert args tuple to a list, so as to have the flexibility or 
	#mutability advantage of modifying the list later on, if needed:
	cpg_file_list = list(cpgs)
	x = pd.read_csv(cpg_file_list[0], sep="\t", names=["chrom","start","end", "perc_meth", "meth", "unmeth"])
	#x.columns = ["chrom","start","end", "meth", "unmeth"]

	for cpg_file in cpg_file_list[1:]:
		current_file = pd.read_csv(cpg_file, sep="\t", names=["chrom","start","end", "perc_meth", "meth", "unmeth"])
		x = pd.merge(x, current_file, **kwargs)
	return x

merged_df = merge_all(cpg_file_1, cpg_file_2, cpg_file_3, how='inner', on=['chrom', 'start', 'end'] )
#merged_df.to_csv(output_file, sep="\t",header=True, index= False)




#line_splitted = open(cleaned_peak_file_1.fn).read().split("\n")
# for each in line_splitted:
#     tab_splitted = each.split()
#     chrom = tab_splitted[0]
#     start = tab_splitted[1]
#     stop = tab_splitted[2]

# count = 0
# for each in line_splitted:
#     if each is not "" or None:
#         print each.split()
#         count +=1
# print count

# c= a_with_b
# d= c.groupby(g=[1, 4, 5], c=10, ops=['sum']
# #c = a_with_b.saveas('/home/surya/surya/for_gata_parse/pybed_overlap.bed', trackline='track name="a and b"')
# c = a_with_b.saveas('/home/surya/surya/for_gata_parse/pybed_overlap.bed')
#print(c)



###
#Cpg_file reading by pandas_dataframe
#####
def load_stringf_for_pybed(input_file):

	df = pd.read_csv(input_file, sep= "\t", header= None)
	df.columns = ["chrom", "start", "end", ".","meth", ".", "unmeth"]
	df.loc[:,"end"] = df["end"] + 1
	df.loc[:,"strand"] = "."
	df.loc[:,"blank_2"] = "."
	Cpg_file = df.loc[:,["chrom","start", "end", "meth", "unmeth","blank_2"]]
	Cpg_bed_filter = Cpg_file.loc[(df["chrom"] != "*")]

	Cpg_list = Cpg_bed_filter.values.tolist()
	Cpg_string_list = [ "\t".join(map(str, line)) for line in Cpg_list]
	Cpg_string_bed_format = "\n".join(Cpg_string_list)

	Cpg_bed_file = pybedtools.BedTool(Cpg_string_bed_format,from_string=True)

	return(Cpg_bed_file)

#Cpg_bed_file = load_stringf_for_pybed(input_file)



#Usage instruction: supply annotated refgene file to the main() function, and
#rest of the output file will be saved to the output dir location pointed out underneath the main()

#Enter the chromosome no. range that you would like to analyse the data on. By default, it would
#take the autosomal gene model from (chr1:chr22, excluding chrX & chrY), however, set chrom_XY=True
#for normal gene model analysis.
#refgene_file = "/home/surya/surya/hudson_alpha_genomics_assignment/refGene_hg19_head.txt"
refgene_file = "/home/surya/refgene_head.txt"
fasta_file = "/home/surya/surya/hudson_alpha_genomics_assignment/hg19-male.fa"

#These would be your output files after you run your script, give the location for output dir:
refgene_model_out_file = "/home/surya/filtered_refgene_model_head.txt"
promoter_coord_file = "/home/surya/promoter_coord_out_head.txt"


def chrom_analysis_range(range_of_chrom_list=22, chrom_XY=False):
	chrom_list = []
	for chrom_no in range(range_of_chrom_list):
		chrom_list.append("chr" + str(chrom_no + 1))
	if chrom_XY:
		chrom_list.append("chrX")
		chrom_list.append("chrY")

	return(chrom_list)



def single_chrom_based_analysis(which_chrom):
	chrom_list = ["chr" + str(which_chrom)]	
	return(chrom_list)



def parse_refgene_model(refgene_file):

	with open(refgene_file, "r") as file:
		sorted_refgene_list = []

		ref_file_read = file.readlines()
		header_splitted = ref_file_read[0].split()
		header_tuple = (header_splitted[2], header_splitted[4], header_splitted[5], header_splitted[3],header_splitted[12])
		header = "\t".join(header_tuple)
		#single_chrom_analysis = single_chrom_based_analysis(19)
		#autosomal_chrom_list = chrom_analysis_range()
		normal_chrom_list = chrom_analysis_range(22, True)

		for line in ref_file_read[1:]:
			splitted = line.strip().split()
			#if str(splitted[2]) in single_chrom_analysis:
			if str(splitted[2]) in normal_chrom_list:
				if splitted[3] == "+":
					sorted_refgene_list.append([str(splitted[2]), int(splitted[4]), int(splitted[5]), str(splitted[3]),str(splitted[12])])
				if splitted[3] == "-":
					sorted_refgene_list.append([str(splitted[2]), int(splitted[5]), int(splitted[4]), str(splitted[3]), str(splitted[12])])

		return(sorted_refgene_list,header)


parsed_refgene = parse_refgene_model(refgene_file)


def final_refgene_model(refgene_file, output_file):

	sorted_refgene_list, header = parse_refgene_model(refgene_file)
	df_refgene = pd.DataFrame(sorted_refgene_list, columns = header.split())
	grouped_refgene = df_refgene.groupby(["chrom","strand","name2"], sort = False, as_index=False)
	refgene_model = grouped_refgene.agg({"txStart": np.min,"txEnd": np.max})	
 	colsorted_final_model = refgene_model.reindex(columns=header.split())
	colsorted_final_model.to_csv(output_file, sep="\t", header = True, index= False)

	return(colsorted_final_model)


final_refgene_model_df = final_refgene_model(refgene_file, refgene_model_out_file)
refgene_list = final_refgene_model_df.values.tolist()



#Creating the dictionary from pandas dataframe for final refgene model to create the bins for each tss in the next step.
#########
refgene_dict_pandas = final_refgene_model_df.to_dict()
header = refgene_dict_pandas.keys()
motif_data = []
for i,name in enumerate(header):
    motif_data.append(refgene_dict_pandas[name].values())

analysis_dict_pandas = dict(zip(header,motif_data))

#In [1081]: analysis_dict_pandas.keys()
#Out[1081]: ['chrom', 'name2', 'txEnd', 'strand', 'txStart']
chrom = analysis_dict_pandas["chrom"]
tss_start = analysis_dict_pandas["txStart"]
tx_end = analysis_dict_pandas["txEnd"]
strand = analysis_dict_pandas["strand"]
gene_name_2 = analysis_dict_pandas["name2"]



### For refseq(tss) analysis using fastinterval and bedtools logic, so this is zero-based so the end element is exclusive:
bin_size = 100
upstream = 1000
downstream = 1000
bins = range(-upstream, (downstream), bin_size)


test_data = [] 
analysis_dict = {}
final_data = list()

#gene_name, chr, strand, start, end, bin_start, bin_end, group
final_header = ["chrom", "chrom_start", "chrom_end", "bin_start", "tss_midpt", "bin_end", "gene_name", "strand"]
print_header = "\t".join(final_header)

for i in range(len(final_header)):
	final_data.append(list())

coordinate_file = "/home/surya/surya/for_gata_parse/tss_coordinate.bed"
with open(coordinate_file, "w") as outfile:
	print print_header
	for i in range(len(chrom)):
		chrome = chrom[i]
		#print chrome
		tss_midpt = tss_start[i]
		strand1 = strand[i]
		gene_name1 = gene_name_2[i]

		for each in bins :
			bin_start = each
			bin_end = (each + bin_size)
	
			if strand1 == "+":
				chrom_start = tss_midpt + (bin_start)
				chrom_end = tss_midpt + (bin_end)
			
			if strand1 == "-":
				#Meth_r
				chrom_start = tss_midpt - (bin_end )
				chrom_end = tss_midpt - (bin_start)

				#This gives the real/original results for the -ve strand, which is true that upstream coordinates 
				#would be greater than the  tss, and the downstream coordinates would be smaller than the tss.
				#But for the sake of bedtools, or other standard formats, its important that we have the stop/end site
				#to be greater than the start site. Thus, it would be obtained by flipping the end site to start site from
				# the following code. However, this could be easily obtained by switching the bin_start side to bin_end
				# as presented in the above code without having to switch to anything. 
				#chrom_start = tss_midpt - (bin_start)
				#chrom_end = tss_midpt - (bin_end)

			line_coords = [chrome, str(chrom_start), str(chrom_end), str(bin_start), str(tss_midpt), str(bin_end), gene_name1, str(strand1)]
			test_data.append("\t".join(line_coords))

			if chrome in analysis_dict:
				analysis_dict[chrome].append((chrome, int(chrom_start), int(chrom_end), gene_name1, "\t".join(line_coords)))
			else:
				analysis_dict[chrome] = [(chrome, int(chrom_start), int(chrom_end), gene_name1, "\t".join(line_coords))]

				
			for i,item in enumerate(line_coords):
				final_data[i].append(item)

			print_coords = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(chrome, chrom_start, chrom_end, bin_start, tss_midpt, bin_end, strand1, gene_name1) 
			print print_coords
			outfile.write(print_coords + "\n")

	final_zip = zip(final_header, final_data)
	final_dict = dict(final_zip)

#####
get_refseq_df_pd = pd.DataFrame(final_dict)
#Uncomment below, if you are looking for one based coordinate overlap like in methplot.R:
#Else, follow normal one that is zero based coordinate overlap like in interval_tree, or bedtools:

#get_refseq_df_pd["chrom_end"] = get_refseq_df_pd["chrom_end"].astype(int) + 1
get_refseq_df = get_refseq_df_pd.loc[:, final_header]
refseq_list = get_refseq_df.values.tolist()

refseq_string_list = [ "\t".join(map(str, line)) for line in refseq_list]
refseq_string_bed_format = "\n".join(refseq_string_list)

refseq_bed_file = pybedtools.BedTool(refseq_string_bed_format,from_string=True)


#Cpg_file reading by pandas_dataframe
#####
df = pd.read_csv("/home/surya/surya/for_gata_parse/SL60583_GATA2_rep1_0_based_1.cov", sep= "\t", header= None)
df.columns = ["chrom", "start", "end", ".","meth", ".", "unmeth"]
df.loc[:,"end"] = df["end"]
df.loc[:,"strand"] = "."
df.loc[:,"blank_2"] = "."
Cpg_file = df.loc[:,["chrom","start", "end", "meth", "unmeth","blank_2"]]
Cpg_bed_filter = Cpg_file.loc[(df["chrom"] != "*")]
# Cpg_bed_filter.to_csv("/home/surya/surya/for_gata_parse/SL60583_GATA2_rep1_pd.cov", sep = "\t", header = False, index = False)
# Cpg_file = "/home/surya/surya/for_gata_parse/SL60583_GATA2_rep1_pd.cov"
# Cpg_bed_file = pybedtools.BedTool(Cpg_file)
# Cpg_with_refseq_test = Cpg_bed_file.intersect(refseq_bed_file, wa = True, wb = True)

Cpg_list = Cpg_bed_filter.values.tolist()
Cpg_string_list = [ "\t".join(map(str, line)) for line in Cpg_list]
Cpg_string_bed_format = "\n".join(Cpg_string_list)

Cpg_bed_file = pybedtools.BedTool(Cpg_string_bed_format,from_string=True)
Cpg_with_refseq = Cpg_bed_file.intersect(refseq_bed_file, wa = True, wb = True)




###Coverting bed_file format/ bedtools output format or any other [[list in list]] format to pandas df:
#######################
Cpg_header = ['chrom', 'start', 'stop','meth', 'unmeth', "strand"]
Cpg_header.extend(final_header)

with open(Cpg_with_refseq.fn) as file:
    bed_list_comp = [ line.strip().split("\t") for line in file ]

Cpg_bed_list_pd = pd.DataFrame(bed_list_comp, columns = Cpg_header)
Cpg_bed_list_pd[["meth", "unmeth"]] =  Cpg_bed_list_pd[["meth", "unmeth"]].astype(int)
Cpg_bed_list_pd.dtypes
#Cpg_bed_list_pd_1[["meth", "unmeth"]] =  Cpg_bed_list_pd_1.loc[:,["meth", "unmeth"]].apply(lambda x: pd.to_numeric(x, errors='coerce'))

#Also, could use this directly, to load the dataframe from bedtool output:
#Cpg_bed_list_pd.columns = Cpg_header
# df = pd.read_table(Cpg_with_refseq.fn, header = Cpg_header)

Cpg_grouped = Cpg_bed_list_pd.groupby(["bin_start", "bin_end"])


directory = "/home/surya/Desktop/scripts/data/groupby_data"
if not os.path.exists(directory):
    os.makedirs(directory)
    #os.rmdir(prior_dir)
    #shutil.rmtree()


#Check if file exists to prevent the append mode for each keys/ch_state value being written to the file:
file_name_list = file_name_list = ["_".join(list(key)) + ".txt" if type(key) == tuple else key + ".txt" for key in Cpg_grouped.groups.keys()]
for file_name in file_name_list:
	prior_file = directory + "/" +file_name
	if os.path.exists(prior_file):
		os.remove(prior_file) 



#print without group_keys:
#for key in Cpg_grouped.groups.keys():
for key in Cpg_grouped.groups:
	grouped_data = Cpg_grouped.get_group(key)

	if type(key) == tuple:
		key = "_".join(list(key))
	else:
		key = key

	#name the output file with its keyword:
	file_name = key + ".txt"
	outfile_path = directory + "/" + file_name
	grouped_data.to_csv(outfile_path, sep="\t", header=True, index = False)
	
	print grouped_data
    #print Cpg_grouped.get_group(key)

#Same as above, print with group_keys:
for each in Cpg_grouped:
	print each

f = {"meth": np.sum, "unmeth": np.sum}
Cpg_agg = Cpg_grouped.agg(f)
Cpg_agg["Percent_meth"] = (Cpg_agg["meth"]) / (Cpg_agg["meth"] + Cpg_agg["unmeth"] ) * 100


Cpg_agg_2 = Cpg_agg.reset_index()
Cpg_agg_2["bin_start"] = Cpg_agg_2["bin_start"].astype(int)
Cpg_agg_2.sort_values(["bin_start"], ascending = [1])
#Cpg_agg.sort_values("bin_start")

Cpg_agg_2.reset_index().to_csv(directory + "/" + "group_data_tss_refseq.txt", sep ="\t", header = True, index = False)




###########################################
###########################################
###########################################
###########################################
###########################################

import numpy as np
import pandas as pd
import pybedtools
import re
import os

#with open("/home/surya/head.txt", "r") as file:

#refgene_file = "/home/surya/refgene_head.txt"
fasta_file = "/home/surya/surya/hudson_alpha_genomics_assignment/hg19-male.fa"

#These would be your output files after you run your script, give the location for output dir:
refgene_model_out_file = "/home/surya/filtered_refgene_model_head.txt"
promoter_coord_file = "/home/surya/promoter_coord_out_head.txt"



### For motif analysis using fastinterval and bedtools logic, so this is zero-based so the end element is exclusive:
bin_size = 100
upstream = 1000
downstream = 1000
bins = range(-upstream, (downstream), bin_size)


test_data = [] 
analysis_dict = {}

#for i in range(5):
#    print analysis_dict["chr1"][i][3]

final_data = list()
final_header = ["chrom", "chrom_start", "chrom_end", "bin_start", "mid_point", "bin_end"]
print_header = "\t".join(final_header)

for i in range(len(final_header)):
	final_data.append(list())

output_file = "/home/surya/surya/for_gata_parse/gata_motif_peak_coords.bed"
input_file = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/Pol2_rep1_peak_chr1_head.bed" 
with open(input_file, "r") as file:
	with open(output_file, "w") as outfile:
		print print_header
		for line in file:
			splitted = line.strip().split()
			chrome = splitted[0]
			#print chrome
			motif_start = int(splitted[1])
			motif_end = int(splitted[2])
			mid_point = (motif_start + motif_end) / 2
			#sequence = sequence_motif_pyfasta[i]

			for each in bins :
				bin_start = each
				bin_end = (each + bin_size)

				chrom_start = mid_point + bin_start
				chrom_end = mid_point + bin_end

				line_coords = [chrome, str(chrom_start), str(chrom_end), str(bin_start), str(mid_point), str(bin_end)]
				test_data.append("\t".join(line_coords))

				if chrome in analysis_dict:
					analysis_dict[chrome].append((chrome, int(chrom_start), int(chrom_end), "\t".join(line_coords)))
				else:
					analysis_dict[chrome] = [(chrome, int(chrom_start), int(chrom_end), "\t".join(line_coords))]

					
				for i,item in enumerate(line_coords):
					final_data[i].append(item)

				print_coords = "%s\t%s\t%s\t%s\t%s\t%s" %(chrome, chrom_start, chrom_end, bin_start, mid_point, bin_end) 
				print print_coords
				outfile.write(print_coords + "\n")

		final_zip = zip(final_header, final_data)
		final_dict = dict(final_zip)


#####
get_refseq_df_pd = pd.DataFrame(final_dict)
#Uncomment below, if you are looking for one based coordinate overlap like in methplot.R:
#Else, follow normal one that is zero based coordinate overlap like in interval_tree, or bedtools:

#get_refseq_df_pd["chrom_end"] = get_refseq_df_pd["chrom_end"].astype(int) + 1
get_refseq_df = get_refseq_df_pd.loc[:, final_header]
refseq_list = get_refseq_df.values.tolist()

refseq_string_list = [ "\t".join(map(str, line)) for line in refseq_list]
refseq_string_bed_format = "\n".join(refseq_string_list)

refseq_bed_file = pybedtools.BedTool(refseq_string_bed_format,from_string=True)



#Cpg file reading by pandas dataframe:
#####
df = pd.read_csv("/home/surya/surya/for_gata_parse/SL60583_GATA2_rep1_0_based_1.cov", sep= "\t", header= None)
df.columns = ["chrom", "start", "end", ".","meth", ".", "unmeth"]
df.loc[:,"end"] = df["end"]
df.loc[:,"strand"] = "."
df.loc[:,"blank_2"] = "."
Cpg_file = df.loc[:,["chrom","start", "end", "meth", "unmeth","blank_2"]]
Cpg_bed_filter = Cpg_file.loc[(df["chrom"] != "*")]
# Cpg_bed_filter.to_csv("/home/surya/surya/for_gata_parse/SL60583_GATA2_rep1_pd.cov", sep = "\t", header = False, index = False)
# Cpg_file = "/home/surya/surya/for_gata_parse/SL60583_GATA2_rep1_pd.cov"
# Cpg_bed_file = pybedtools.BedTool(Cpg_file)
# Cpg_with_refseq_test = Cpg_bed_file.intersect(refseq_bed_file, wa = True, wb = True)

Cpg_list = Cpg_bed_filter.values.tolist()
Cpg_string_list = [ "\t".join(map(str, line)) for line in Cpg_list]
Cpg_string_bed_format = "\n".join(Cpg_string_list)

Cpg_bed_file = pybedtools.BedTool(Cpg_string_bed_format,from_string=True)
Cpg_with_refseq = Cpg_bed_file.intersect(refseq_bed_file, wa = True, wb = True)




###Coverting bed_file format/ bedtools output format or any other [[list in list]] format to pandas df:
#######################
Cpg_header = ['chrom', 'start', 'stop','meth', 'unmeth', "strand"]
Cpg_header.extend(final_header)

with open(Cpg_with_refseq.fn) as file:
    bed_list_comp = [ line.strip().split("\t") for line in file ]

Cpg_bed_list_pd = pd.DataFrame(bed_list_comp, columns = Cpg_header)
Cpg_bed_list_pd[["meth", "unmeth"]] =  Cpg_bed_list_pd[["meth", "unmeth"]].astype(int)
Cpg_bed_list_pd.dtypes
#Cpg_bed_list_pd_1[["meth", "unmeth"]] =  Cpg_bed_list_pd_1.loc[:,["meth", "unmeth"]].apply(lambda x: pd.to_numeric(x, errors='coerce'))

#Also, could use this directly, to load the dataframe from bedtool output:
#Cpg_bed_list_pd.columns = Cpg_header
# df = pd.read_table(Cpg_with_refseq.fn, header = Cpg_header)

Cpg_grouped = Cpg_bed_list_pd.groupby(["bin_start", "bin_end"])


directory = "/home/surya/Desktop/scripts/data/groupby_data"
if not os.path.exists(directory):
    os.makedirs(directory)
    #os.rmdir(prior_dir)
    #shutil.rmtree()


#Check if file exists to prevent the append mode for each keys/ch_state value being written to the file:
file_name_list = file_name_list = ["_".join(list(key)) + ".txt" if type(key) == tuple else key + ".txt" for key in Cpg_grouped.groups.keys()]
for file_name in file_name_list:
	prior_file = directory + "/" +file_name
	if os.path.exists(prior_file):
		os.remove(prior_file) 



#print without group_keys:
#for key in Cpg_grouped.groups.keys():
for key in Cpg_grouped.groups:
	grouped_data = Cpg_grouped.get_group(key)

	if type(key) == tuple:
		key = "_".join(list(key))
	else:
		key = key

	#name the output file with its keyword:
	file_name = key + ".txt"
	outfile_path = directory + "/" + file_name
	grouped_data.to_csv(outfile_path, sep="\t", header=True, index = False)
	
	print grouped_data
    #print Cpg_grouped.get_group(key)

#Same as above, print with group_keys:
for each in Cpg_grouped:
	print each

f = {"meth": np.sum, "unmeth": np.sum}
Cpg_agg = Cpg_grouped.agg(f)
Cpg_agg["Percent_meth"] = (Cpg_agg["meth"]) / (Cpg_agg["meth"] + Cpg_agg["unmeth"] ) * 100


Cpg_agg_2 = Cpg_agg.reset_index()
Cpg_agg_2["bin_start"] = Cpg_agg_2["bin_start"].astype(int)
Cpg_agg_2.sort_values(["bin_start"], ascending = [1])
#Cpg_agg.sort_values("bin_start")

Cpg_agg_2.reset_index().to_csv(directory + "/" + "group_data__bed_peak.txt", sep ="\t", header = True, index = False)





cpg_file_1 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/k562_wgbs_chr1.cov"
cpg_file_2 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/SL60583_GATA2_rep1.cov_chr1.cov"
cpg_file_3 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/SL63653_GATA2_rep2.cov_chr1.cov"


output_dir = "/home/surya/Desktop/scripts/data/groupby_data_test"
ggplot_outfile_suffix = "group_data_motifs_concat.txt"



#def capture_binned_meth_perc_list(*meth_file, **kwargs):


meth_file = [cpg_file_1, cpg_file_2, cpg_file_3]
meth_file_list = list(meth_file)

frames = []

for idx, each_file in enumerate(meth_file_list):
	cpg_bed_file = load_stringf_for_cpg_after_edits(each_file)
	meth_data = generate_motifs_binned_perc_meth(motif_bed_file, cpg_bed_file, output_dir, ggplot_outfile_suffix)
	meth_data["file_name"] = os.path.splitext(basename(each_file))[0]
	frames.append(meth_data)

appended_data = pd.concat(frames, ignore_index=True)



##############################
#############################
#############################
#############################




def generate_motifs_binned_perc_meth(motif_bedfile, cpg_bedfile, **kwargs): 
	
	print "kwargs: ", kwargs

	Cpg_with_motif = cpg_bedfile.intersect(motif_bedfile, wa = True, wb = True)	
	final_header = ["chrom", "chrom_start", "chrom_end", "bin_start", "mid_point", "bin_end"]

	###Coverting bed_file format/ bedtools output format or any other [[list in list]] format to pandas df:
	Cpg_header = ['chrom', 'start', 'stop','meth', 'unmeth', "strand"]
	Cpg_header.extend(final_header)
	with open(Cpg_with_motif.fn) as file:
	    bed_list_comp = [ line.strip().split("\t") for line in file ]
	Cpg_bed_list_pd = pd.DataFrame(bed_list_comp, columns = Cpg_header)
	Cpg_bed_list_pd[["meth", "unmeth"]] =  Cpg_bed_list_pd[["meth", "unmeth"]].astype(int)
	Cpg_grouped = Cpg_bed_list_pd.groupby(["bin_start", "bin_end"])

	directory = kwargs["output_dir"]
	if not os.path.exists(directory):
	    os.makedirs(directory)

		#Check if file exists to prevent the append mode for each keys/ch_state value being written to the file:
	file_name_list = [kwargs["file_basename"] + "_" + "_".join(list(key)) + ".txt" if type(key) == tuple else kwargs["file_basename"] + "_" + key + ".txt" for key in Cpg_grouped.groups.keys()]
	for file_name in file_name_list:
		prior_file = directory + "/" +file_name
		if os.path.exists(prior_file):
			os.remove(prior_file) 

	for key in Cpg_grouped.groups:
		grouped_data = Cpg_grouped.get_group(key)

		if type(key) == tuple:
			key = "_".join(list(key))
		else:
			key = key

		#name the output file with its keyword:
		file_name = kwargs["file_basename"] + "_" + key + ".txt"
		outfile_path = directory + "/" + file_name
		print outfile_path
		grouped_data.to_csv(outfile_path, sep="\t", header=True, index = False)	


	f = {"meth": np.sum, "unmeth": np.sum}
	Cpg_agg = Cpg_grouped.agg(f)
	Cpg_agg["Percent_meth"] = (Cpg_agg["meth"]) / (Cpg_agg["meth"] + Cpg_agg["unmeth"] )
	#Cpg_agg["group_name"] = ggplot_group_name

	Cpg_agg_1 = Cpg_agg.reset_index()
	Cpg_agg_1["bin_start"] = Cpg_agg_1["bin_start"].astype(int)
	Cpg_agg_sorted = Cpg_agg_1.sort_values(["bin_start"], ascending = [1])
	#Cpg_agg_sorted.index = Cpg_agg_sorted.index.sort_values()

	Cpg_agg_sorted.to_csv(directory + "/" + kwargs["ggplot_outfile"], sep ="\t", header = True, index = False)

	return(Cpg_agg_sorted)




cpg_file_1 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/k562_wgbs_chr1.cov"
cpg_file_2 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/SL60583_GATA2_rep1.cov_chr1.cov"
cpg_file_3 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/SL63653_GATA2_rep2.cov_chr1.cov"


def capture_binned_meth_perc_list(meth_file_list, **kwargs):

	
	frames = []
	for idx, each_file in enumerate(meth_file_list):
		file_basename = os.path.splitext(basename(each_file))[0]	
		cpg_bed_file = load_stringf_for_cpg_after_edits(each_file)
		meth_data = generate_motifs_binned_perc_meth(motif_bed_file, cpg_bed_file, file_basename=file_basename, **kwargs)
		meth_data["file_name"] = os.path.splitext(basename(each_file))[0]
		frames.append(meth_data)
	appended_data = pd.concat(frames, ignore_index=True)
	#appended_data.to_csv()

	return(appended_data)

list_of_files = [cpg_file_2]
plot_data_set = capture_binned_meth_perc_list(list_of_files, output_dir=output_dir, ggplot_outfile=ggplot_outfile_suffix)

In [41]: import glob

In [42]: a = glob.glob("/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/*chr1.cov")

In [43]: a
Out[43]: 
['/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/MACS_SL60583_SL60585_10_peaks.bed_chr1.cov',
 '/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/MACS_SL63653_SL63655_10_peaks.bed_chr1.cov',
 '/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/wgbs_GATA2.cov_chr1.cov',
 '/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/SL60583_GATA2_rep1.cov_chr1.cov',
 '/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/SL63653_GATA2_rep2.cov_chr1.cov']


appended_data = []
for infile in glob.glob("*.xlsx"):
    data = pandas.read_excel(infile)
    appended_data.append(data) ## store dataframes in list
appended_data = pd.concat(appended_data, axis=1) # equivalent to join operation

appended_data .to_excel('appedned.xlsx')



#meth_file = [cpg_file_1, cpg_file_2, cpg_file_3]

# def capture_binned_meth_perc_list(*meth_file, **kwargs):

# 	meth_file_list = list(meth_file)
# 	frames = []

# 	for idx, each_file in enumerate(meth_file_list):
# 		file_basename = os.path.splitext(basename(each_file))[0]	
# 		cpg_bed_file = load_stringf_for_cpg_after_edits(each_file)
# 		meth_data = generate_motifs_binned_perc_meth(motif_bed_file, cpg_bed_file, file_basename=file_basename, **kwargs)
# 		meth_data["file_name"] = os.path.splitext(basename(each_file))[0]
# 		frames.append(meth_data)
# 	appended_data = pd.concat(frames, ignore_index=True)
# 	#appended_data.to_csv()

# 	return(appended_data)

# plot_data_set = capture_binned_meth_perc_list(cpg_file_2, cpg_file_3, output_dir=output_dir, ggplot_outfile=ggplot_outfile_suffix, multiple_files=True)





# cpg_bed_file = load_stringf_for_cpg_after_edits(meth_file_list[0])
# meth_data_1 = generate_motifs_binned_perc_meth(motif_bed_file, cpg_bed_file, output_dir, ggplot_outfile_suffix, "k562_gata_rep1")

# for each_file in meth_file_list[1:]:
# 	cpg_bed_file = load_stringf_for_cpg_after_edits(each_file)
# 	meth_data = generate_motifs_binned_perc_meth(motif_bed_file, cpg_bed_file, output_dir, ggplot_outfile_suffix, "k562_gata_rep1")
# 	meth_data["file_name"] = os.path.splitext(os.path.basename(each_file))[0]
# 	meth_data_1.append(meth_data)









#Multi-indexing with index name as "bin_start" and "bin_end":
In [186]: Cpg_agg = Cpg_grouped.agg({"meth": np.sum, "unmeth": np.sum})
In [276]: Cpg_agg["Percent_meth"] = (Cpg_agg["meth"]) / (Cpg_agg["meth"] + Cpg_agg["unmeth"] ) * 100

In [277]: Cpg_agg
Out[277]: In [81]: tss_binned_perc_meth
Out[81]: 
    bin_start bin_end  meth  unmeth  Percent_meth
0        -100       0     0     105      0.000000
1       -1000    -900    18      66     21.428571
2        -200    -100     1     104      0.952381
3        -300    -200    10      78     11.363636
4        -400    -300     5      62      7.462687
5        -500    -400    21      61     25.609756
6        -600    -500    12      70     14.634146
7        -700    -600    15      47     24.193548
8        -800    -700     5      22     18.518519
9        -900    -800     5      10     33.333333
10          0     100     1     145      0.684932
11        100     200     4     120      3.225806
12        200     300     0      83      0.000000
13        300     400     6      22     21.428571
14        400     500     1      43      2.272727
15        500     600     1      36      2.702703
16        600     700     3      31      8.823529
17        700     800     7      42     14.285714
18        800     900     5      42     10.638298
19        900    1000     6      39     13.333333


                   meth  unmeth  Percent_meth
bin_start bin_end                            
-100      0           0     105      0.000000
-1000     -900       18      66     21.428571
-200      -100        1     104      0.952381
-300      -200       10      79     11.235955
-400      -300        5      62      7.462687
-500      -400       21      61     25.609756
-600      -500       12      70     14.634146
-700      -600       15      47     24.193548
-800      -700        5      22     18.518519
-900      -800        5      10     33.333333
0         100         1     145      0.684932
100       200         4     120      3.225806
200       300         0      83      0.000000
300       400         6      24     20.000000
400       500         1      43      2.272727
500       600         1      36      2.702703
600       700         3      36      7.692308
700       800         7      43     14.000000
800       900         5      42     10.638298
900       1000        7      43     14.000000

In [213]: Cpg_agg.index
Out[213]: 
MultiIndex(levels=[[u'-100', u'-1000', u'-200', u'-300', u'-400', u'-500', u'-600', u'-700', u'-800', u'-900', u'0', u'100', u'200', u'300', u'400', u'500', u'600', u'700', u'800', u'900'], [u'-100', u'-200', u'-300', u'-400', u'-500', u'-600', u'-700', u'-800', u'-900', u'0', u'100', u'1000', u'200', u'300', u'400', u'500', u'600', u'700', u'800', u'900']],
           labels=[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], [9, 8, 0, 1, 2, 3, 4, 5, 6, 7, 10, 12, 13, 14, 15, 16, 17, 18, 19, 11]],
           names=[u'bin_start', u'bin_end'])


In [214]: Cpg_agg.keys()
Out[214]: Index([u'meth', u'unmeth'], dtype='object')


In [190]: Cpg_agg.columns
Out[190]: Index([u'meth', u'unmeth'], dtype='object')


In [252]: Cpg_agg.index.values
Out[252]: 
array([('-100', '0'), ('-1000', '-900'), ('-200', '-100'),
       ('-300', '-200'), ('-400', '-300'), ('-500', '-400'),
       ('-600', '-500'), ('-700', '-600'), ('-800', '-700'),
       ('-900', '-800'), ('0', '100'), ('100', '200'), ('200', '300'),
       ('300', '400'), ('400', '500'), ('500', '600'), ('600', '700'),
       ('700', '800'), ('800', '900'), ('900', '1000')], dtype=object)

In [250]: list(Cpg_agg.index.values)
Out[250]: 
[('-100', '0'),
 ('-1000', '-900'),
 ('-200', '-100'),
 ('-300', '-200'),
 ('-400', '-300'),
 ('-500', '-400'),
 ('-600', '-500'),
 ('-700', '-600'),
 ('-800', '-700'),
 ('-900', '-800'),
 ('0', '100'),
 ('100', '200'),
 ('200', '300'),
 ('300', '400'),
 ('400', '500'),
 ('500', '600'),
 ('600', '700'),
 ('700', '800'),
 ('800', '900'),
 ('900', '1000')]



########

#Single index or indexing with the index name as "bin_start":
In [274]: Cpg_agg.reset_index("bin_start")
Out[274]: 
        bin_start  meth  unmeth  Percent_meth
bin_end                                      
0            -100     0     105      0.000000
-900        -1000    18      66     21.428571
-100         -200     1     104      0.952381
-200         -300    10      79     11.235955
-300         -400     5      62      7.462687
-400         -500    21      61     25.609756
-500         -600    12      70     14.634146
-600         -700    15      47     24.193548
-700         -800     5      22     18.518519
-800         -900     5      10     33.333333
100             0     1     145      0.684932
200           100     4     120      3.225806
300           200     0      83      0.000000
400           300     6      24     20.000000
500           400     1      43      2.272727
600           500     1      36      2.702703
700           600     3      36      7.692308
800           700     7      43     14.000000
900           800     5      42     10.638298
1000          900     7      43     14.000000



In [211]: Cpg_agg.reset_index("bin_start").index
Out[211]: 
Index([u'0', u'-900', u'-100', u'-200', u'-300', u'-400', u'-500', u'-600',
       u'-700', u'-800', u'100', u'200', u'300', u'400', u'500', u'600',
       u'700', u'800', u'900', u'1000'],
      dtype='object', name=u'bin_end')


In [208]: Cpg_agg.reset_index("bin_start").columns
Out[208]: Index([u'bin_start', u'meth', u'unmeth'], dtype='object')

In [210]: Cpg_agg.reset_index("bin_start").keys()
Out[210]: Index([u'bin_start', u'meth', u'unmeth'], dtype='object')



#######

#Normal range indexing, since all the index name has been dissolved, however, range indexing is the most fundamental indexing,
#so, it would always remain there, with no changes.
#Basically, to make it short, changes to default range indexing like in normal dataframes:

In [273]: Cpg_agg.reset_index()
Out[273]: 
   bin_start bin_end  meth  unmeth  Percent_meth
0       -100       0     0     105      0.000000
1      -1000    -900    18      66     21.428571
2       -200    -100     1     104      0.952381
3       -300    -200    10      79     11.235955
4       -400    -300     5      62      7.462687
5       -500    -400    21      61     25.609756
6       -600    -500    12      70     14.634146
7       -700    -600    15      47     24.193548
8       -800    -700     5      22     18.518519
9       -900    -800     5      10     33.333333
10         0     100     1     145      0.684932
11       100     200     4     120      3.225806
12       200     300     0      83      0.000000
13       300     400     6      24     20.000000
14       400     500     1      43      2.272727
15       500     600     1      36      2.702703
16       600     700     3      36      7.692308
17       700     800     7      43     14.000000
18       800     900     5      42     10.638298
19       900    1000     7      43     14.000000


In [204]: Cpg_agg.reset_index().index
Out[204]: RangeIndex(start=0, stop=20, step=1)

In [212]: Cpg_agg.reset_index().columns
Out[212]: Index([u'bin_start', u'bin_end', u'meth', u'unmeth'], dtype='object')

In [205]: Cpg_agg.reset_index().keys()
Out[205]: Index([u'bin_start', u'bin_end', u'meth', u'unmeth'], dtype='object')



#Set the index to the normal dataframe:
In [217]: a = Cpg_bed_list_pd.head()

In [218]: a
Out[218]: 
  chrom    start     stop  meth  unmeth strand chrom chrom_start chrom_end  \
0  chr1  1050735  1050737     0       3      .  chr1     1050736   1050836   
1  chr1  1050737  1050739     0       3      .  chr1     1050736   1050836   
2  chr1  1050755  1050757     0       1      .  chr1     1050736   1050836   
3  chr1  1050756  1050758     0       1      .  chr1     1050736   1050836   
4  chr1  1050761  1050763     0       1      .  chr1     1050736   1050836   

  bin_start tss_midpt bin_end gene_name strand  
0       900   1051736    1000  C1orf159      -  
1       900   1051736    1000  C1orf159      -  
2       900   1051736    1000  C1orf159      -  
3       900   1051736    1000  C1orf159      -  
4       900   1051736    1000  C1orf159      -  

In [219]: a.index
Out[219]: RangeIndex(start=0, stop=5, step=1)


In [220]: a.set_index("chrom")
Out[220]: 
                start     stop  meth  unmeth strand chrom_start chrom_end  \
chrom                                                                       
(chr1, chr1)  1050735  1050737     0       3      .     1050736   1050836   
(chr1, chr1)  1050737  1050739     0       3      .     1050736   1050836   
(chr1, chr1)  1050755  1050757     0       1      .     1050736   1050836   
(chr1, chr1)  1050756  1050758     0       1      .     1050736   1050836   
(chr1, chr1)  1050761  1050763     0       1      .     1050736   1050836   

             bin_start tss_midpt bin_end gene_name strand  
chrom                                                      
(chr1, chr1)       900   1051736    1000  C1orf159      -  
(chr1, chr1)       900   1051736    1000  C1orf159      -  
(chr1, chr1)       900   1051736    1000  C1orf159      -  
(chr1, chr1)       900   1051736    1000  C1orf159      -  
(chr1, chr1)       900   1051736    1000  C1orf159      -  

In [221]: a.set_index("chrom").index
Out[221]: 
Index([(u'chr1', u'chr1'), (u'chr1', u'chr1'), (u'chr1', u'chr1'),
       (u'chr1', u'chr1'), (u'chr1', u'chr1')],
      dtype='object', name=u'chrom')

In [222]: a.set_index("chrom")["start"]
Out[222]: 
chrom
(chr1, chr1)    1050735
(chr1, chr1)    1050737
(chr1, chr1)    1050755
(chr1, chr1)    1050756
(chr1, chr1)    1050761
Name: start, dtype: object



In [244]: Cpg_agg.reset_index("bin_start")
Out[244]: 
        bin_start  meth  unmeth  Percent_meth
bin_end                                      
0            -100     0     105      0.000000
-900        -1000    18      66     21.428571
-100         -200     1     104      0.952381
-200         -300    10      79     11.235955
-300         -400     5      62      7.462687
-400         -500    21      61     25.609756
-500         -600    12      70     14.634146
-600         -700    15      47     24.193548
-700         -800     5      22     18.518519
-800         -900     5      10     33.333333
100             0     1     145      0.684932
200           100     4     120      3.225806
300           200     0      83      0.000000
400           300     6      24     20.000000
500           400     1      43      2.272727
600           500     1      36      2.702703
700           600     3      36      7.692308
800           700     7      43     14.000000
900           800     5      42     10.638298
1000          900     7      43     14.000000



# Changing the index value to the list from numpy to python list:
In [245]: Cpg_agg.reset_index("bin_start").index
Out[245]: 
Index([u'0', u'-900', u'-100', u'-200', u'-300', u'-400', u'-500', u'-600',
       u'-700', u'-800', u'100', u'200', u'300', u'400', u'500', u'600',
       u'700', u'800', u'900', u'1000'],
      dtype='object', name=u'bin_end')

In [246]: Cpg_agg.reset_index("bin_start").index.values
Out[246]: 
array(['0', '-900', '-100', '-200', '-300', '-400', '-500', '-600', '-700',
       '-800', '100', '200', '300', '400', '500', '600', '700', '800',
       '900', '1000'], dtype=object)

In [247]: list(Cpg_agg.reset_index("bin_start").index.values)
Out[247]: 
['0',
 '-900',
 '-100',
 '-200',
 '-300',
 '-400',
 '-500',
 '-600',
 '-700',
 '-800',
 '100',
 '200',
 '300',
 '400',
 '500',
 '600',
 '700',
 '800',
 '900',
 '1000']














###
#awk '{print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$6}' SL60583_GATA2_rep1.cov > SL60583_GATA2_rep1_0_based.cov
import pybedtools

#Cpg_file = "/home/surya/surya/for_gata_parse/head_cpg.bed"
Cpg_file = "/home/surya/surya/for_gata_parse/SL60583_GATA2_rep1_0_based_1.cov"
coordinate_file = "/home/surya/surya/for_gata_parse/gata_motif.bed"
c = "/home/surya/surya/for_gata_parse/fastinterval_overlaps.bed"

cpg_bed = pybedtools.example_bedtool(Cpg_file)
gata_bed = pybedtools.example_bedtool(coordinate_file)
c = pybedtools.example_bedtool(c)
a_with_b = cpg_bed.intersect(gata_bed, wa = True, wb = True)
#a_with_b = a.intersect(a, u = True)
#a_with_b = cpg_bed.intersect(gata_bed, wo= True)
#a_with_b.groupby(g=[1, 4, 5], c=10, o =[sum])
a_with_b_grouped = a_with_b.groupby(g=[11, 15, 16], c=[5,7])
c = a_with_b_grouped.saveas('/home/surya/surya/for_gata_parse/grouped_gata_motif.bed')


a_with_b.groupby(g=[11, 15, 16], c=[5,7]).head()
a_with_b.head()

a_with_b.groupby(g=[11, 15, 16], c=[5,7]).tail()
a_with_b.tail()
 #df = pandas.read_table(cpg_bed.fn, names=['chrom', 'start', 'stop', 'perc', 'meth', 'strand', 'unmeth'])
#print len(a_with_b)
#print(a_with_b)
#print(open(a_with_b.fn).read())
d= groupby(g=[1, 4, 5], c=[5,10], o=["sum","max"])

cpg_intersect = a.intersect( 
								[b.fn, c.fn],
								names = ["b", "c"],
								filenames = True,
								wa = True,
								u = True

							)
cpg_intersect.head()



#line_splitted = open(a_with_b.fn).read().split("\n")
# for each in line_splitted:
#     tab_splitted = each.split()
#     chrom = tab_splitted[0]
#     start = tab_splitted[1]
#     stop = tab_splitted[2]

count = 0
for each in line_splitted:
    if each is not "" or None:
        print each.split()
        count +=1
print count

c= a_with_b
d= groupby(g=[1, 4, 5], c=10, ops=['sum'])
#print d

#c = a_with_b.saveas('/home/surya/surya/for_gata_parse/pybed_overlap.bed', trackline='track name="a and b"')
c = a_with_b_grouped.saveas('/home/surya/surya/for_gata_parse/grouped_gata_motif.bed')



######################


def extract_tss_promoter_coordinates(final_model_df, bp_upstream, bp_downstream, output_file):
	header = ["chrom","txStart","strand","name2"]
	get_tss_df = final_model_df.loc[:,header]
	get_tss_df["Start"] = 0
	get_tss_df["End"] = 0

	plus_strand = (get_tss_df.strand == "+")
	get_tss_df["Start"][plus_strand] = get_tss_df.txStart + ( - int(bp_upstream))
	get_tss_df["End"][plus_strand] = get_tss_df.txStart + ( + int(bp_downstream))

	minus_strand = (get_tss_df.strand == "-")
	get_tss_df["Start"][minus_strand] = get_tss_df.txStart - ( + int(bp_downstream))
	get_tss_df["End"][minus_strand] = get_tss_df.txStart - ( - int(bp_upstream) )
	#Eliminate start coordinate < 0:
	get_tss_df = get_tss_df[get_tss_df.Start > 0]
	select_cols = ["chrom","Start","End","strand","name2"]
	colsorted_tss_df = get_tss_df.reindex(columns=select_cols)
	get_tss_python_list = colsorted_tss_df.values.tolist()
	colsorted_tss_df.to_csv(output_file, sep="\t",header=False, index= False)

	return(colsorted_tss_df, get_tss_python_list)



get_tss_df, get_tss_python_list = extract_tss_promoter_coordinates(final_refgene_model_df, 1000, 50, promoter_coord_file)


########################
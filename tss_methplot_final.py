import numpy as np
import pandas as pd
import pybedtools
import re
import os
import time
from os.path import basename

start_time = time.time()


input_file_1 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/pol2/Pol2_rep1_peak_chr1.bed")
input_file_2 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/pol2/Pol2_rep2_peak_chr1.bed")
input_file_3 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/gata/MACS_SL60583_SL60585_10_peaks.bed")
input_file_4 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/gata/MACS_SL63653_SL63655_10_peaks.bed")

#refgene_file = os.path.expanduser("/home/surya/surya/hudson_alpha_genomics_assignment/refGene_hg19_head.txt"
fasta_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/hg_19_misc/hg19-male.fa")
refgene_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/hg_19_misc/refGene_head.txt")
cpg_file_1 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/k562_wgbs_chr1.cov")
cpg_file_2 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/gata/SL60583_GATA2_rep1.cov_chr1.cov")
cpg_file_3 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/gata/SL63653_GATA2_rep2.cov_chr1.cov")


output_dir = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/groupby_refseq_final")
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

ggplot_outfile_suffix = "group_data_motifs_concating_tss.txt"
tss_coordinate_out_file = output_dir + "/" + "tss_coordinate_out.bed"
refgene_model_out_file = output_dir + "/" + "filtered_refgene_model_head.txt"



def intersect_all(*beds, **kwargs):

	#convert args tuple to a list, so as to have the flexibility, or mutability advantage of modifying the list later on, if needed:
	bedlist = list(beds)
	x = pybedtools.BedTool(bedlist[0])

	for bed in bedlist[1:]:
		x = x.intersect(bed, **kwargs)
	return x


#f(x):cleaned_peak_file_1 = intersect_all(input_file_1, input_file_2, input_file_3, input_file_4, wa = True, u = True)
#cleaned_peak_file_1.saveas('/home/surya/surya/for_gata_parse/pybed_overlap.bed')



#Merging all the cpg files with the common key (chrom, start, end):
def merge_all(*cpgs, **kwargs):

	#convert args tuple to a list, so as to have the flexibility, or mutability advantage of modifying the list later on, if needed:
	cpg_file_list = list(cpgs)
	x = pd.read_csv(cpg_file_list[0], sep="\t", names=["chrom","start","end", "perc_meth", "meth", "unmeth"])

	for cpg_file in cpg_file_list[1:]:
		current_file = pd.read_csv(cpg_file, sep="\t", names=["chrom","start","end", "perc_meth", "meth", "unmeth"])
		x = pd.merge(x, current_file, **kwargs)
	return x


#f(x):merged_df = merge_all(cpg_file_1, cpg_file_2, cpg_file_3, how='inner', on=['chrom', 'start', 'end'] )
#merged_df.to_csv(output_file, sep="\t",header=True, index= False)



#Enter the chromosome no. range that you would like to analyse the data on. By default, it would take the autosomal 
#gene model from (chr1:chr22, excluding chrX & chrY), however, set chrom_XY=True for normal gene model analysis.
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



def final_refgene_model(refgene_file, output_file):

	sorted_refgene_list, header = parse_refgene_model(refgene_file)
	df_refgene = pd.DataFrame(sorted_refgene_list, columns = header.split())
	grouped_refgene = df_refgene.groupby(["chrom","strand","name2"], sort = False, as_index=False)
	refgene_model = grouped_refgene.agg({"txStart": np.min,"txEnd": np.max})	
 	colsorted_final_model = refgene_model.reindex(columns=header.split())
	colsorted_final_model.to_csv(output_file, sep="\t", header = True, index= False)

	return(colsorted_final_model)



def generate_dict_from_pandas_df_dict(final_refgene_model_df):

	#Creating the dictionary from pandas dataframe for final refgene model to create the bins for each tss in the next step.
	#refgene_list = final_refgene_model_df.values.tolist()
	refgene_dict_pandas = final_refgene_model_df.to_dict()
	header = refgene_dict_pandas.keys()
	motif_data = []

	for i,name in enumerate(header):
	    motif_data.append(refgene_dict_pandas[name].values())

	analysis_dict_pandas = dict(zip(header,motif_data))
	chrom_1 = analysis_dict_pandas["chrom"]
	tss_start_1 = analysis_dict_pandas["txStart"]
	tx_end_1 = analysis_dict_pandas["txEnd"]
	strand_1 = analysis_dict_pandas["strand"]
	gene_name_1 = analysis_dict_pandas["name2"]

	return(chrom_1, tss_start_1, tx_end_1, strand_1, gene_name_1)

#coordinates_info_for_binning = generate_dict_from_pandas_df_dict(final_refgene_model_df_pd) 
#chrom, tss_start, tx_end, strand, gene_name_2 = generate_dict_from_pandas_df_dict(final_refgene_model_df_pd)


def generate_tss_binned_coords(coordinates_info, upstream_range, downstream_range, bin_size, coord_outfile):
	### For refseq(tss) analysis using fastinterval and bedtools logic, so this is zero-based so the end element is exclusive:
	binSize = bin_size
	upstream = upstream_range
	downstream = downstream_range
	bins = range(-upstream, (downstream), binSize)

	test_data = [] 
	analysis_dict = {}
	final_data = list()

	#gene_name, chr, strand, start, end, bin_start, bin_end, group
	final_header = ["chrom", "chrom_start", "chrom_end", "bin_start", "tss_midpt", "bin_end", "gene_name", "strand"]
	print_header = "\t".join(final_header)

	for i in range(len(final_header)):
		final_data.append(list())

	with open(coord_outfile, "w") as outfile:
		#print print_header

		for i in range(len(coordinates_info[0])):
			chrome = coordinates_info[0][i]
			#print chrome
			tss_midpt = coordinates_info[1][i]
			strand1 = coordinates_info[3][i]
			gene_name1 = coordinates_info[4][i]

			for each in bins :
				bin_start = each
				bin_end = (each + bin_size)
		
				if strand1 == "+":
					chrom_start = tss_midpt + (bin_start)
					chrom_end = tss_midpt + (bin_end)
				
				if strand1 == "-":
					#Meth_r concept, start and end switched, so as to always maintain the higher coords for end site
					chrom_start = tss_midpt - (bin_end )
					chrom_end = tss_midpt - (bin_start)

					#for original data w/o switching the start to end for -ve strand, 
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
				#print print_coords
				outfile.write(print_coords + "\n")

		final_zip = zip(final_header, final_data)
		final_dict = dict(final_zip)

	return(final_dict)




def load_stringf_for_refseq_from_flydict(final_dict):

	final_header = ["chrom", "chrom_start", "chrom_end", "bin_start", "tss_midpt", "bin_end", "gene_name", "strand"]
	get_refseq_df_pd = pd.DataFrame(final_dict)
	
	#Uncomment below, if you are looking for one based coordinate overlap like in methplot.R, else, follow normal one that is zero based coordinate overlap like in interval_tree, or bedtools:
	#get_refseq_df_pd["chrom_end"] = get_refseq_df_pd["chrom_end"].astype(int) + 1
	
	get_refseq_df = get_refseq_df_pd.loc[:, final_header]
	refseq_list = get_refseq_df.values.tolist()

	refseq_string_list = [ "\t".join(map(str, line)) for line in refseq_list]
	refseq_string_bed_format = "\n".join(refseq_string_list)
	refseq_bed_file = pybedtools.BedTool(refseq_string_bed_format,from_string=True)

	return(refseq_bed_file)




#Cpg_file reading via pandas_dataframe for adding the strands for standard pybed usage, and editing the coordinates:
def load_stringf_for_cpg_after_edits(input_file):

	df = pd.read_csv(input_file, sep= "\t", header= None)
	df.columns = ["chrom", "start", "end", "perc_meth", "meth", "unmeth"]
	df.loc[:,"end"] = df["end"] + 1
	df.loc[:,"strand"] = "."
	Cpg_file = df.loc[:,["chrom","start", "end", "meth", "unmeth", "strand"]]
	#Cpg_bed_filter = Cpg_file.loc[(df["chrom"] != "*")]

	#Cpg_list = Cpg_bed_filter.values.tolist()
	Cpg_list = Cpg_file.values.tolist()
	Cpg_string_list = [ "\t".join(map(str, line)) for line in Cpg_list]
	Cpg_string_bed_format = "\n".join(Cpg_string_list)

	Cpg_bed_file = pybedtools.BedTool(Cpg_string_bed_format,from_string=True)

	return(Cpg_bed_file)


#f(x):cpg_bed_file = load_stringf_for_cpg_after_edits(cpg_input_file)


#Cpg_file reading via pandas_dataframe for adding the strands for standard pybed usage, and editing the coordinates:
def load_pybedObject_for_cpg_after_edits(input_file):

	file_basename = os.path.splitext(basename(input_file))[0]
	file_ext = os.path.splitext(basename(input_file))[1]

	df = pd.read_csv(input_file, sep= "\t", header= None)
	df.columns = ["chrom", "start", "end", "perc_meth", "meth", "unmeth"]
	df.loc[:,"end"] = df["end"] + 1
	df.loc[:,"strand"] = "."
	Cpg_file = df.loc[:,["chrom","start", "end", "meth", "unmeth", "strand"]]
	#Cpg_bed_filter = Cpg_file.loc[(df["chrom"] != "*")]
	
	file_name = output_dir + "/" + file_basename + "_" + "cpg_edit" + file_ext
	#Make sure that header=False else, pybedtools won't recognise the headers above:
	Cpg_file.to_csv(file_name, sep ="\t", header = False, index = False)
	Cpg_bed_file = pybedtools.BedTool(file_name)

	#Cpg_list = Cpg_bed_filter.values.tolist()
	# Cpg_list = Cpg_file.values.tolist()
	# Cpg_string_list = [ "\t".join(map(str, line)) for line in Cpg_list]
	# Cpg_string_bed_format = "\n".join(Cpg_string_list)
	#Cpg_bed_file = pybedtools.BedTool(Cpg_string_bed_format,from_string=True)

	return(Cpg_bed_file)


#f(x):cpg_bed_file = load_pybedObject_for_cpg_after_edits(cpg_input_file)



def generate_tss_binned_perc_meth(refseq_bedfile, cpg_bedfile, **kwargs):

	print "kwargs: ", kwargs		
	Cpg_with_refseq = cpg_bedfile.intersect(refseq_bedfile, wa = True, wb = True)	
	final_header = ["chrom", "chrom_start", "chrom_end", "bin_start", "tss_midpt", "bin_end", "gene_name", "strand"]

	###Coverting bed_file format/ bedtools output format or any other [[list in list]] format to pandas df:
	Cpg_header = ['chrom', 'start', 'stop','meth', 'unmeth', "strand"]
	Cpg_header.extend(final_header)
	with open(Cpg_with_refseq.fn) as file:
	    bed_list_comp = [ line.strip().split("\t") for line in file ]

	Cpg_bed_list_pd = pd.DataFrame(bed_list_comp, columns = Cpg_header)
	# Cpg_bed_list_grouped = Cpg_ bed_list.groupby(["features"])
	# for each_features in Cpg_bed_list_grouped.groups():
	# 	do the following grouping meth, unmeth......

	Cpg_bed_list_pd[["meth", "unmeth"]] =  Cpg_bed_list_pd[["meth", "unmeth"]].astype(int)
	Cpg_grouped = Cpg_bed_list_pd.groupby(["bin_start", "bin_end"])

	directory = output_dir
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
		
		grouped_data.to_csv(outfile_path, sep="\t", header=True, index = False)	


	f = {"meth": np.sum, "unmeth": np.sum}
	Cpg_agg = Cpg_grouped.agg(f)
	Cpg_agg["Percent_meth"] = (Cpg_agg["meth"]) / (Cpg_agg["meth"] + Cpg_agg["unmeth"] )
	#Cpg_agg["group_name"] = ggplot_group_name

	Cpg_agg_1 = Cpg_agg.reset_index()
	Cpg_agg_1["bin_start"] = Cpg_agg_1["bin_start"].astype(int)
	Cpg_agg_sorted = Cpg_agg_1.sort_values(["bin_start"], ascending = [1])
	#Cpg_agg_sorted.index = Cpg_agg_sorted.index.sort_values()

	Cpg_agg_sorted.to_csv(directory + "/" + ggplot_outfile_suffix, sep ="\t", header = True, index = False)

	return(Cpg_agg_sorted)



#f(x):tss_binned_perc_meth = generate_tss_binned_perc_meth(refseq_bed_file, cpg_bed_file, output_dir=output_dir, ggplot_outfile_suffix=ggplot_outfile_suffix, file_basename="")
#tss_binned_perc_meth.index = tss_binned_perc_meth.index.sort_values()


def capture_binned_meth_perc_list(meth_file_list, refseq_tss_bed_file):	
	frames = []
	for idx, each_file in enumerate(meth_file_list):
		file_basename1 = os.path.splitext(basename(each_file))[0]	
		#cpg_bed_file = load_stringf_for_cpg_after_edits(each_file)
		cpg_bed_file = load_pybedObject_for_cpg_after_edits(each_file)
		meth_data = generate_tss_binned_perc_meth(refseq_tss_bed_file, cpg_bed_file, file_basename=file_basename1)
		meth_data["file_name"] = os.path.splitext(basename(each_file))[0]
		frames.append(meth_data)
	appended_data = pd.concat(frames, ignore_index=True)
	appended_data.to_csv(output_dir + "/" + ggplot_outfile_suffix, sep ="\t", header = True, index = False)

	return(appended_data)


def main():

	#global refseq_bed_file

	#f(x):cleaned_peak_file_1 = intersect_all(input_file_1, input_file_2, input_file_3, input_file_4, wa = True, u = True)
	#cleaned_peak_file_1.saveas('/home/surya/surya/for_gata_parse/pybed_overlap.bed')

	#f(x):merged_df = merge_all(cpg_file_1, cpg_file_2, cpg_file_3, how='inner', on=['chrom', 'start', 'end'] )
	#merged_df.to_csv(output_file, sep="\t",header=True, index= False)


	final_refgene_model_df_pd = final_refgene_model(refgene_file, refgene_model_out_file)
	coordinates_info_for_binning = generate_dict_from_pandas_df_dict(final_refgene_model_df_pd) 

	final_dict_tss_coord = generate_tss_binned_coords(coordinates_info_for_binning, 1000, 1000, 100, tss_coordinate_out_file )
	refseq_bed_file = load_stringf_for_refseq_from_flydict(final_dict_tss_coord)

	list_of_files = [cpg_file_2, cpg_file_3, cpg_file_1]
	plot_data_set = capture_binned_meth_perc_list(list_of_files, refseq_bed_file)



#refseq_bed_file = None

if __name__ == '__main__':
	main()

else:
	print "Functions Imported from other module"


print "Time for analysis = ", time.time()-start_time
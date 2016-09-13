
import numpy as np
import pandas as pd
import pybedtools
import re
import os
from os.path import basename

input_file_1 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/Pol2_rep1_peak_chr1.bed"
input_file_2 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/Pol2_rep2_peak_chr1.bed"
input_file_3 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/MACS_SL60583_SL60585_10_peaks.bed_chr1.cov"
input_file_4 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/MACS_SL63653_SL63655_10_peaks.bed_chr1.cov"


def intersect_all(*beds, **kwargs):

	#convert args tuple to a list, so as to have the flexibility, or mutability advantage of modifying the list later on, if needed:
	bedlist = list(beds)
	x = pybedtools.BedTool(bedlist[0])

	for bed in bedlist[1:]:
		x = x.intersect(bed, **kwargs)
	return x

#f(x):cleaned_peak_file_1 = intersect_all(input_file_1, input_file_2, input_file_3, input_file_4, wa = True, u = True)
#cleaned_peak_file_1.saveas('/home/surya/surya/for_gata_parse/pybed_overlap.bed')



cpg_file_1 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/k562_wgbs_chr1.cov"
cpg_file_2 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/SL60583_GATA2_rep1.cov_chr1.cov"
cpg_file_3 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/SL63653_GATA2_rep2.cov_chr1.cov"
header=["chrom","start","end", "perc_meth", "meth", "unmeth"]


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



peak_input_file = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/Pol2_rep1_peak_chr1_head.bed" 
cpg_file_1 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/k562_wgbs_chr1.cov"
cpg_file_2 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/SL60583_GATA2_rep1.cov_chr1.cov"
cpg_file_3 = "/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/SL63653_GATA2_rep2.cov_chr1.cov"

output_dir = "/home/surya/Desktop/scripts/data/groupby_data_final_peak"
if not os.path.exists(output_dir):
	os.makedirs(output_dir)
	  
ggplot_outfile_suffix = "group_data_motifs_concating.txt"
peak_coord_outfile = output_dir + "/" + "gata_motif_peak_coords.bed"


def generate_motifs_binned_coords(input_file, upstream_range, downstream_range, bin_size, coord_outfile):
	### For motif analysis using fastinterval and bedtools logic, so this is zero-based so the end element is exclusive:
	binSize = bin_size
	upstream = upstream_range
	downstream = downstream_range
	bins = range(-upstream, (downstream), binSize)

	test_data = [] 
	analysis_dict = {}

	final_data = list()
	final_header = ["chrom", "chrom_start", "chrom_end", "bin_start", "mid_point", "bin_end"]
	print_header = "\t".join(final_header)

	for i in range(len(final_header)):
		final_data.append(list())

	with open(input_file, "r") as file:
		with open(coord_outfile, "w") as outfile:
			print print_header
			for line in file:
				splitted = line.strip().split()
				chrome = splitted[0]
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

	return(final_dict)


final_dict_motif_coord = generate_motifs_binned_coords(peak_input_file, 1000, 1000, 100, peak_coord_outfile)



def load_stringf_for_motif_from_flydict(final_dict):

	final_header = ["chrom", "chrom_start", "chrom_end", "bin_start", "mid_point", "bin_end"]
	get_peaks_df_pd = pd.DataFrame(final_dict)	
	#Uncomment below, if you are looking for one based coordinate overlap like in methplot.R, else, follow normal one that is zero based coordinate overlap like in interval_tree, or bedtools:
	#get_peaks_df_pd["chrom_end"] = get_peaks_df_pd["chrom_end"].astype(int) + 1

	get_peaks_df = get_peaks_df_pd.loc[:, final_header]
	peaks_list = get_peaks_df.values.tolist()

	peaks_string_list = [ "\t".join(map(str, line)) for line in peaks_list]
	peaks_string_bed_format = "\n".join(peaks_string_list)
	peaks_bed_file = pybedtools.BedTool(peaks_string_bed_format,from_string=True)

	return(peaks_bed_file)


motif_bed_file = load_stringf_for_motif_from_flydict(final_dict_motif_coord)


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

	df = pd.read_csv(input_file, sep= "\t", header= None)
	df.columns = ["chrom", "start", "end", "perc_meth", "meth", "unmeth"]
	df.loc[:,"end"] = df["end"] + 1
	df.loc[:,"strand"] = "."
	Cpg_file = df.loc[:,["chrom","start", "end", "meth", "unmeth", "strand"]]
	#Cpg_bed_filter = Cpg_file.loc[(df["chrom"] != "*")]
	
	file_name = output_dir + "/" + file_basename + "_" + "cpg_edit"
	#Make sure that header=False else, pybedtools won't recognise the headers above:
	Cpg_file.to_csv(file_name, sep ="\t", header = False, index = False)
	Cpg_bed_file = pybedtools.BedTool(file_name)


	#Cpg_list = Cpg_bed_filter.values.tolist()
	# Cpg_list = Cpg_file.values.tolist()
	# Cpg_string_list = [ "\t".join(map(str, line)) for line in Cpg_list]
	# Cpg_string_bed_format = "\n".join(Cpg_string_list)
	# Cpg_bed_file = pybedtools.BedTool(Cpg_string_bed_format,from_string=True)

	return(Cpg_bed_file)


#f(x):cpg_bed_file = load_pybedObject_for_cpg_after_edits(cpg_input_file)


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
	#Cpg_agg_sorted.to_csv(directory + "/" + ggplot_outfile_suffix, sep ="\t", header = True, index = False)

	return(Cpg_agg_sorted)


def capture_binned_meth_perc_list(meth_file_list):	
	frames = []
	for idx, each_file in enumerate(meth_file_list):
		file_basename = os.path.splitext(basename(each_file))[0]	
		cpg_bed_file = load_stringf_for_cpg_after_edits(each_file)
		meth_data = generate_motifs_binned_perc_meth(motif_bed_file, cpg_bed_file, file_basename=file_basename)
		meth_data["file_name"] = os.path.splitext(basename(each_file))[0]
		frames.append(meth_data)
	appended_data = pd.concat(frames, ignore_index=True)
	appended_data.to_csv(output_dir + "/" + ggplot_outfile_suffix, sep ="\t", header = True, index = False)

	return(appended_data)

list_of_files = [cpg_file_2, cpg_file_3]
plot_data_set = capture_binned_meth_perc_list(list_of_files)


##################
##################
##################
##################


def capture_binned_meth_perc_list(meth_file_list):	
	frames = []
	for idx, each_file in enumerate(meth_file_list):
		file_basename = os.path.splitext(basename(each_file))[0]	
		cpg_bed_file = load_stringf_for_cpg_after_edits(each_file)
		meth_data = generate_motifs_binned_perc_meth(motif_bed_file, cpg_bed_file, file_basename=file_basename)
		meth_data["file_name"] = os.path.splitext(basename(each_file))[0]
		frames.append(meth_data)
	appended_data = pd.concat(frames, ignore_index=True)
	appended_data.to_csv(output_dir + "/" + ggplot_outfile_suffix, sep ="\t", header = True, index = False)

	return(appended_data)


def main():

	global motif_bed_file

	#f(x):cleaned_peak_file_1 = intersect_all(input_file_1, input_file_2, input_file_3, input_file_4, wa = True, u = True)
	#cleaned_peak_file_1.saveas('/home/surya/surya/for_gata_parse/pybed_overlap.bed')

	#f(x):merged_df = merge_all(cpg_file_1, cpg_file_2, cpg_file_3, how='inner', on=['chrom', 'start', 'end'] )
	#merged_df.to_csv(output_file, sep="\t",header=True, index= False)

	final_dict_motif_coord = generate_motifs_binned_coords(peak_input_file, 1000, 1000, 100, peak_coord_outfile)
	motif_bed_file = load_stringf_for_motif_from_flydict(final_dict_motif_coord)

	list_of_cpg_files = [cpg_file_2, cpg_file_3]
	plot_data_set = capture_binned_meth_perc_list(list_of_cpg_files)



motif_bed_file = None

if __name__ == '__main__':
	main()

else:
	print "Functions Imported from other module"


# In [17]: kwargs = {"file_name": "test.txt"}

# In [18]: kwargs
# Out[18]: {'file_name': 'test.txt'}

# In [19]: print "kwargs", kwargs
# kwargs {'file_name': 'test.txt'}

# In [20]: print "kwargs:", kwargs
# kwargs: {'file_name': 'test.txt'}




print "Time for analysis = ", time.time()-start_time
#list_of_files = [path1, path2,...]
	
#I can see: file_list = glob.glob('/path to directory/*[pattern]/*.txt") working for me. 
#Also, sub_dir = [dir[:1] for dir in glob.glob("/path to directory/*[pattern]/*[pattern]/")]
#This is bascially using the Unix shell rule.   
#@User123: that doesn't list directories recursively. You are listing all text files one level deep, but not in further subdirectories or even directly in path to directory.
#Example:
#In [118]: glob.glob("/home/surya/Desktop/scripts/data/*/*.cov")
#Out[118]: ['/home/surya/Desktop/scripts/data/meth_chip_bs_data/k562_wgbs_chr1.cov']

# Here, though meth_chip_bs_data has further subdirectories like gata, and pol2 containing the *.cov file, but this doesn't return all the .cov files
# but only lists all .cov files one level deep, or immediate subdirectory only, and not even the sub directory of the immediate subdirectory. Here eg: gata and pol2
# '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/',

# In [121]: glob.glob("/home/surya/Desktop/scripts/data/*/*/*.cov")
# Out[121]: 
# ['/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/Pol2_rep1_cpg_chr1.cov',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/SL88935_Pol2_rep2.cov',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/wgbs_Pol2.cov',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/Pol2_rep2_cpg_chr1.cov',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/SL68321_Pol2_rep1.cov',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/MACS_SL60583_SL60585_10_peaks.bed_chr1.cov',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/wgbs_GATA2.cov',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/MACS_SL63653_SL63655_10_peaks.bed_chr1.cov',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/wgbs_GATA2.cov_chr1.cov',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/SL60583_GATA2_rep1.cov_chr1.cov',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/SL60583_GATA2_rep1.cov',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/SL63653_GATA2_rep2.cov_chr1.cov',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/SL63653_GATA2_rep2.cov']




# In [61]: a = glob.glob("/home/surya/Desktop/scripts/data/meth_chip_bs_data/*")

# In [62]: a
# Out[62]: 
# ['/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/k562_wgbs_chr1.cov',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/script.sh']

# In [72]: a = [file for file in glob.glob("/home/surya/Desktop/scripts/data/meth_chip_bs_data/*")]

# In [73]: a
# Out[73]: 
# ['/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/k562_wgbs_chr1.cov',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/script.sh']



# open a series of subfolders in a folder and find some text files and print some lines of the text files:
# Current versions of glob.glob() cannot list files in subdirectories recursively (an update that will be included in python 3.5 when released adds a ** option for arbitrary nested directory traversal).

#use os.walk() combined with fnmatch.filter() instead:

# import os
# import fnmatch

# path = 'C:/Users/sam/Desktop/file1'

# configfiles = [os.path.join(dirpath, f)
#     for dirpath, dirnames, files in os.walk(path)
#     for f in fnmatch.filter(files, '*.txt')]

# This'll walk your directories recursively and return all absolute pathnames to matching .txt files. In this specific case the fnmatch.filter() may be overkill, you could also use a .endswith() test:

# import os

# path = 'C:/Users/sam/Desktop/file1'

# configfiles = [os.path.join(dirpath, f)
#     for dirpath, dirnames, files in os.walk(path)
#     for f in files if f.endswith('.txt')]



# In [72]: a = [file for file in glob.glob("/home/surya/Desktop/scripts/data/meth_chip_bs_data/*")]

# In [73]: a
# Out[73]: 
# ['/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/k562_wgbs_chr1.cov',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/script.sh']

# In [74]: a = glob.iglob("/home/surya/Desktop/scripts/data/*2/*")

# In [75]: a
# Out[75]: <generator object iglob at 0x7f6085fd2690>

# In [76]: a = glob.glob("/home/surya/Desktop/scripts/data/*2/*")

# In [77]: a
# Out[77]: []

# In [78]: a = glob.glob("/home/surya/Desktop/scripts/data/meth_chip_bs_data/*2/*")

# In [79]: a
# Out[79]: 
# ['/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/Pol2_rep1_cpg_chr1.cov',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/SL88935_Pol2_rep2.cov',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/MACS_SL88935_SL88941_15_peaks.bed',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/Pol2_rep1_peak_chr1_head.bed',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/wgbs_Pol2.cov',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/Pol2_rep2_cpg_chr1.cov',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/Pol2_rep2_peak_chr1.bed',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/SL68321_Pol2_rep1.cov',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/Pol2_rep1_peak_chr1.bed',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/MACS_SL68321_SL68325_15_peaks.bed']

# In [80]: a = glob.glob("/home/surya/Desktop/scripts/data/meth_chip_bs_data/*2/*.bed")

# In [81]: a
# Out[81]: 
# ['/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/MACS_SL88935_SL88941_15_peaks.bed',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/Pol2_rep1_peak_chr1_head.bed',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/Pol2_rep2_peak_chr1.bed',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/Pol2_rep1_peak_chr1.bed',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/MACS_SL68321_SL68325_15_peaks.bed']

# In [82]: a = glob.glob("/home/surya/Desktop/scripts/data/meth_chip_bs_data/*ta/*.bed")

# In [83]: a
# Out[83]: 
# ['/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/MACS_SL60583_SL60585_10_peaks.bed',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/MACS_SL63653_SL63655_10_peaks.bed']

# In [84]: a = glob.glob("/home/surya/Desktop/scripts/data/meth_chip_bs_data/*/*.bed")

# In [85]: a
# Out[85]: 
# ['/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/MACS_SL88935_SL88941_15_peaks.bed',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/Pol2_rep1_peak_chr1_head.bed',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/Pol2_rep2_peak_chr1.bed',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/Pol2_rep1_peak_chr1.bed',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/pol2/MACS_SL68321_SL68325_15_peaks.bed',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/MACS_SL60583_SL60585_10_peaks.bed',
#  '/home/surya/Desktop/scripts/data/meth_chip_bs_data/gata/MACS_SL63653_SL63655_10_peaks.bed']



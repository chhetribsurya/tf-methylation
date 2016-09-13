import numpy
import pandas
import pybedtools
import re

#with open("/home/surya/head.txt", "r") as file:


# In [32]: for line in tss_list:
#     test1 = "\t".join(map(str, line))

#Usage instruction: supply annotated refgene file and hg19 fasta file to the main() function, and
#rest of the output file will be saved to the output dir location pointed out underneath the main()



#Enter the chromosome no. range that you would like to analyse the data on. By default, it would
#take the autosomal gene model from (chr1:chr22, excluding chrX & chrY), however, set chrom_XY=True
#for normal gene model analysis.

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
	df_refgene = pandas.DataFrame(sorted_refgene_list, columns = header.split())
	grouped_refgene = df_refgene.groupby(["chrom","strand","name2"], sort = False, as_index=False)
	refgene_model = grouped_refgene.agg({"txStart": numpy.min,"txEnd": numpy.max})	
 	colsorted_final_model = refgene_model.reindex(columns=header.split())
	colsorted_final_model.to_csv(output_file, sep="\t", header = True, index= False)

	return(colsorted_final_model)




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


# header = pd.DataFrame(dict_1 ).to_dict().keys()
# motif_data = []
# for i,name in enumerate(header):
#     motif_data.append(pd.DataFrame(dict_1 ).to_dict()[name].values())
# dict(zip(header,motif_data))



def extract_possible_TATAbox_coordinates(final_model_df, bp_upstream_1, bp_upstream_2, output_file):
	header = ["chrom","txStart","strand","name2"]
	get_tatabox_df = final_model_df.loc[:,header]
	get_tatabox_df["Start"] = 0
	get_tatabox_df["End"] = 0

	plus_strand = (get_tatabox_df.strand == "+")
	get_tatabox_df["Start"][plus_strand] = get_tatabox_df.txStart + ( - int(bp_upstream_2))
	get_tatabox_df["End"][plus_strand] = get_tatabox_df.txStart + ( - int(bp_upstream_1))

	minus_strand = (get_tatabox_df.strand == "-")
	get_tatabox_df["Start"][minus_strand] = get_tatabox_df.txStart - ( - int(bp_upstream_1))
	get_tatabox_df["End"][minus_strand] = get_tatabox_df.txStart - ( - int(bp_upstream_2) )
	#Eliminate start coordinate < 0:
	get_tatabox_df = get_tatabox_df[get_tatabox_df.Start > 0]
	select_cols = ["chrom","Start","End","strand","name2"]
	colsorted_tatabox_df = get_tatabox_df.reindex(columns=select_cols)
	get_tatabox_python_list = colsorted_tatabox_df.values.tolist()
	colsorted_tatabox_df.to_csv(output_file, sep="\t",header=False, index= False)

	return(colsorted_tatabox_df, get_tatabox_python_list)

# In [378]: grouped_df = df.groupby(["Sepal.Width","Petal.Length"], as_index = False)
# df_a = grouped_df.get_group(key)
# df.is_copy = False (after the indexing operation)
#df_a[df_a.index ==11]
#df_a.index => [11,22,59]
#df.loc[[11,22,59]]
#or, df_a.loc[df_a.index]
#df.iloc[[0,1,2]]
# .loc[row_indexer,col_indexer] = value instead #To troubleshoot warnings:

# df_a.loc[:,"count"] = df_a.loc[:,"Sepal.Length"] + df_a.loc[:,"Petal.Length"]

# In [369]: df_a["count"] = df_a.loc[:,"Sepal.Length"] + df_a.loc[:,"Petal.Length"]

# In [370]: df_a["counts"] = df_a.loc[:,"Sepal.Length"] + df_a.loc[:,"Petal.Length"]

# In [371]: df_a["counting"]  = df_a["Sepal.Length"]  + df_a["Petal.Length"]
# cols = ['col1', 'col2', 'col3']
# df[cols] = df[cols].applymap(some_function)

# result = df\
#     .loc[:, df.columns!='Time']\
#     .groupby(lambda x: x[0], axis=1)\
#     .apply(lambda x: np.sqrt((x**2).sum(1)))\
#     .apply(lambda x: x / df['Time'])

# print result

#             A          B
# 1    1.404626   1.310639
# 2   -2.954644 -10.874091
# 3    3.479836   6.105961
# 4    3.885530   2.244544
# 5    0.995012   1.434228
# 6   11.278208  11.454466
# 7   -1.209242  -1.281165
# 8   -5.175911  -5.905070
# 9   11.889318  16.758958
# 10  -0.978014  -0.590767

# A possible start.

# Filtering out column names corresponding to a particular vector. For example

# In [20]: filter(lambda x: x.startswith("A_"),df.columns)
# Out[20]: ['A_x', 'A_y', 'A_z']
# Sub selecting these columns from the DataFrame

# In [22]: df[filter(lambda x: x.startswith("A_"),df.columns)]
# Out[22]: 
#          A_x       A_y       A_z
# 1  -0.123527 -0.547239 -0.453707
# 2  -0.112098 -1.122609  0.218538
# 3   0.818276 -1.563931  0.097377
# 4   0.766515 -0.650482 -0.087203
# 5  -2.419204 -0.882383  0.005204
# 6   1.115630  0.081825 -1.038442
# 7   0.226386  0.039879  0.732611
# 8   0.389148  0.158289  0.440282
# 9  -0.308314 -0.839347 -0.517989
# 10  0.473552  0.059428  0.726088
# So, using this technique you can get chunks of 3 columns. For example.

# column_initials = ["A","B"]
# for column_initial in column_initials:
#     df["Velocity_"+column_initial]=df[filter(lambda x: x.startswith(column_initial+"_"),df.columns)].apply(lambda x: np.sqrt(x.dot(x)), axis=1)/df.Time


# In [32]: df[['Velocity_A','Velocity_B']]
# Out[32]: 
#     Velocity_A  Velocity_B
# 1    -9.555311  -22.467965
# 2    -5.568487   -7.177625
# 3    -9.086257  -12.030091
# 4     2.007230    1.144208
# 5     1.824531    0.775006
# 6     1.472305    2.623467
# 7     1.954044    3.967796
# 8    -0.485576   -1.384815
# 9    -7.736036   -6.722931
# 10    1.392823    5.369757

# These are techniques to apply function to element, column or dataframe.

# Map: It iterates over each element of a series. 
# df[‘column1’].map(lambda x: 10+x), this will add 10 to each element of column1.
# df[‘column2’].map(lambda x: ‘AV’+x), this will concatenate “AV“ at the beginning of each element of column2 (column format is string).

# Apply: As the name suggests, applies a function along any axis of the DataFrame.
# df[[‘column1’,’column2’]].apply(sum), it will returns the sum of all the values of column1 and column2.

# ApplyMap: This helps to apply a function to each element of dataframe.
# func = lambda x: x+2
# df.applymap(func), it will add 2 to each element of dataframe (all columns of dataframe must be numeric type)


# In [443]: grouped_df = df.groupby(["Sepal.Width","Petal.Length"], as_index= False)
# In [443]: df_1 = grouped_df.agg({"Sepal.Length" : np.sum})
# In [443]: df_1["key"] = 0
# In [443]: df_1["Perc_meth"] = (df_1["Sepal.Length"] + df_1["Petal.Length"])/100

grouped_df = df.groupby(["Sepal.Width","Petal.Length"], as_index= False)
grouped_df = df.groupby(["Sepal.Width","Petal.Length"]).count()
df.groupby(["Sepal.Width","Petal.Length"], as_index= False).apply(funct)
grouped_df.apply(lambda x : x["Sepal.Length"] + 2, axis = 1)
grouped_df.agg(lambda x: x.sum())

df_1 = grouped_df.agg({"Sepal.Length" : np.sum, "Petal.Width" : np.sum})

In [277]: np.random.random((2,3)).tolist()
In [276]: for i in range(5):
    print np.random.randint(1,10)



def parse_hg19_genome_fasta(hg19_fasta_file):
	with open(hg19_fasta_file,"r") as file:
		fasta_file = file.read()
		hg19_fasta_dict = {}

		for each_bed in fasta_file.split(">"):
		    if each_bed and each_bed != "":
		        chrom,seq = each_bed.strip("\n").split("\n")
		        if chrom not in hg19_fasta_dict:
                		hg19_fasta_dict[chrom] = [(seq[1:50])]
                    else:
                        hg19_fasta_dict[chrom].append((seq[1:50]))

    	return(hg_19_fasta_dict) 

#hg19_fasta_dictionary = parse_hg19_genome_fasta(fasta_file)

def extract_seq_and_count_CpG_in_promoter(promoter_coordinates_file, hg_19_fasta_file, output_file):
	with open(output_file, "w") as out_file:
		read_promoter_coord = pybedtools.BedTool(promoter_coordinates_file)
		fasta = pybedtools.example_filename(hg_19_fasta_file)
		extract_fasta = read_promoter_coord.sequence(fi=fasta)
		prom_bedtofasta = open(read_promoter_coord.seqfn).read()
		bed_dict = {}

		print "\n\n\n\nHere is the head of file exhibiting CpG count in promoter region\n..........."
		file_header = "chrom\tstart\tend\tCpG_counts\n"
		print file_header
		out_file.write(file_header)

		for each_bed in prom_bedtofasta.split(">"):
			if each_bed and each_bed != "":
		                meta, seq = each_bed.split("\n", 1)
		                chrom, start, end = re.split("[:-]", meta) 
		                count_CpG= seq.count("CG")
		                save_output = "%s\t%s\t%s\t%s\n" %(chrom, start, end, count_CpG)
		                out_file.write(save_output)
		                print save_output
						                   
		                #Group chromosome wise:
		                if chrom not in bed_dict:
		                    bed_dict[chrom] = [(start,end, count_CpG)]
		                else:
		                    bed_dict[chrom].append((start,end,count_CpG))

		return(bed_dict)




def extract_seq_and_count_TATA_box_in_promoter(tatabox_coord_file, hg_19_fasta_file, output_file):
	with open(output_file, "w") as out_file:
		read_tatabox_coord = pybedtools.BedTool(tatabox_coord_file)
		fasta = pybedtools.example_filename(hg_19_fasta_file)
		extract_fasta = read_tatabox_coord.sequence(fi=fasta)
		tatabox_bedtofasta = open(read_tatabox_coord.seqfn).read()
		bed_dict = {}

		print "\n\n\nHere is the head of file exhibiting TATAbox count in promoter region\n..........."
		file_header = "chrom\tstart\tend\tTATAbox_counts\n"
		print file_header
		out_file.write(file_header)

		for each_bed in tatabox_bedtofasta.split(">"):
			if each_bed and each_bed != "":
		                meta, seq = each_bed.split("\n", 1)
		                chrom, start, end = re.split("[:-]", meta) 
		                count_CpG= seq.count("TATAA")
		                save_output = "%s\t%s\t%s\t%s\n" %(chrom, start, end, count_CpG)
		                out_file.write(save_output)
		                print save_output
						                   
		                #Group chromosome wise:
		                if chrom not in bed_dict:
		                    bed_dict[chrom] = [(start,end, count_CpG)]
		                else:
		                    bed_dict[chrom].append((start,end,count_CpG))

		return(bed_dict)



def combine_CpG_TATAbox_info_file(promoter_CpG_info_file, promoter_TATA_info_file, output_file):
	promoter_CpG_df = pandas.read_csv(promoter_CpG_info_file, sep="\t", header = True, names = "chrom\tstart\tend\tCpG_counts".split())
	promoter_TATA_df = pandas.read_csv(promoter_TATA_info_file, sep="\t", header = True, names = "chrom\tstart\tend\tTATAbox_counts".split())
	combined_CpG_tata_info = pandas.merge(promoter_CpG_df, promoter_TATA_df, how='inner', on=['chrom', 'start', 'end'])
	combined_CpG_tata_info.to_csv(output_file, sep="\t",header=True, index= False)

	return(combined_CpG_tata_info)



def main():	
	# Enter the location for refgene and fasta file:
	refgene_file = "/home/surya/surya/hudson_alpha_genomics_assignment/refGene_hg19"
	fasta_file = "/home/surya/surya/hudson_alpha_genomics_assignment/hg19-male.fa"

	#These would be your output files after you run your script, give the location for output dir:
	refgene_model_out_file = "/home/surya/surya/hudson_alpha_genomics_assignment/test/filtered_refgene_model.txt"
	promoter_coord_file = "/home/surya/surya/hudson_alpha_genomics_assignment/test/promoter_coord_out.txt"
	tatabox_coord_file = "/home/surya/surya/hudson_alpha_genomics_assignment/test/tatabox_coord_out.txt"
	CpG_count_file = "/home/surya/surya/hudson_alpha_genomics_assignment/test/prom_region_CpG_count.txt"
	TATAbox_count_file = "/home/surya/surya/hudson_alpha_genomics_assignment/test/prom_region_TATA_count.txt"
	promoter_CpG_TATAbox_file = "/home/surya/surya/hudson_alpha_genomics_assignment/test/promoter_CpG_TATA3_count.txt"

	#These are the function that would be called after the main() function runs:
	parsed_refgene = parse_refgene_model(refgene_file)
	final_refgene_model_df = final_refgene_model(refgene_file, refgene_model_out_file)
	get_tss_df, get_tss_python_list = extract_tss_promoter_coordinates(final_refgene_model_df, 1000, 50, promoter_coord_file)
	get_tatabox_df, get_tatabox_python_list = extract_possible_TATAbox_coordinates(final_refgene_model_df, 10, 40 , tatabox_coord_file)
	CpG_info_containing_file = extract_seq_and_count_CpG_in_promoter(promoter_coord_file, fasta_file, CpG_count_file)
	TATA_info_containing_file = extract_seq_and_count_TATA_box_in_promoter(promoter_coord_file, fasta_file, TATAbox_count_file)
	combined_CpG_tatabox_file = combine_CpG_TATAbox_info_file(CpG_count_file, TATAbox_count_file, promoter_CpG_TATAbox_file)
	print "\nScript run completed\n\n"



if __name__ == '__main__':
	main()

else:
	print "Functions Imported from other module"




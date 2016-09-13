
import re

#with open("/home/surya/head.txt", "r") as file: 

meme_file = "/home/surya/Desktop/MeMe/meme.txt"
query_motif = "Motif 2"


##########
##########


# Returns the list of motif_regex/motif_seq
def find_motif_regex(MeMe_file, motif_name=False):

    with open(MeMe_file, "r") as motif_file:
        data= motif_file.read()

        motif_name_pattern = "Motif\s\d+\sregular\sexpression"
        regex = re.compile(motif_name_pattern)
        motif_name_list = regex.findall(data)

        motif_seq_pattern = "Motif\s\d+\sregular\sexpression\n-*\n(.*?\n)"
        regex_next = re.compile(motif_seq_pattern)
        motif_seq_list = regex_next.findall(data)

        motif_dict = {}

        for i in range(len(motif_name_list)):
            motif_key, motif_value = " ".join(motif_name_list[i].split()[0:2]), motif_seq_list[i].strip("\n")
            
            #If you don't want list style dictionary, just use the down below commented one,
            #i.e ignoring if and else condition to facilitate the formation of string based value of key straight away.
            #motif_dict[motif_key] = motif_value

            if motif_key in motif_dict:
                motif_dict[motif_key].append(motif_seq_list)
            else:
                motif_dict[motif_key] = [motif_value]

    	#return(motif_dict[query_motif])
    	return(motif_dict)

motif_pattern = find_motif_regex(meme_file)
gata_regex = motif_pattern[query_motif]


########
########



content = open("/home/surya/Desktop/MeMe/meme.txt", "r").read()
#pattern  = r"Motif 2 sites sorted by position p-value.*\n.*\n.*\n.*\n"
pattern  = r"Motif \d+ sites sorted by position p-value.*\n.*\n.*\n.*\n"
regex = re.compile(pattern)
matched_list = regex.findall(content)
header = ["chrom","start","end","strand","sequence"]
motif_header = "\t".join(header)

#$$$
Master_dict = dict()
#$$$

# print matched_list
for each in matched_list:

	#$$$
	M_key = each[0:8].strip()
	print M_key
	#$$$

	motif_data = list()	
	for i in range(len(header)):
		motif_data.append(list())



	#print " -------------------********************** ----------------- ", each[:8]
	el = len(each)
	start_index = content.find(each) + el
	end_index = content.find("-----", start_index)
	required_each_content = content[start_index: end_index]
	each_motif_lines = required_each_content.strip().split("\n")
	
	for motif_line in each_motif_lines:
		splitted = motif_line.split()
		chrom1 = splitted[0].split(":")[0]
		start1 = int(splitted[0].split(":")[1]) + int(splitted[2])
		strand1 = splitted[1]
		sequence1 = splitted[5]
		end1 = start1 + len(sequence1)


		req = [chrom1, start1, end1, strand1, sequence1]

		for i,item in enumerate(req):
			motif_data[i].append(item)

	motif_zip = zip(header,motif_data)
	motif_dict = dict(motif_zip)

	#$$$
	Master_dict[M_key] = motif_dict
	#$$$


test = Master_dict["Motif 2"]
test.keys()
chrom = test["chrom"]
start = test["start"]
end = test["end"]
strand = test["strand"]
sequence_motif = test["sequence"]

##Joining the dictionary:
###################

df_1 = pd.DataFrame(Master_dict["Motif 10"])
df_2 = pd.DataFrame(Master_dict["Motif 18"])
df_3 = pd.DataFrame(Master_dict["Motif 19"])
df_concat = pd.concat([df_1, df_2, df_3])


###################



dict_1 = Master_dict["Motif 10"]
dict_2 = Master_dict["Motif 18"]
dict_3 = Master_dict["Motif 19"]

dict_items = [dict_1, dict_2, dict_3]

final_dict = {}

for each_dict in dict_items:
    for key, value in each_dict.iteritems():
        if key in final_dict:
            final_dict[key].extend(value)
        else:
            final_dict[key] = []
            final_dict[key].extend(value)

df_test = pd.DataFrame(final_dict)




###################
print("Final data starts here...\n\n")
count = 0
for i in range(len(chrom)):
	motif_coordinates = "%s\t%s\t%s\t%s\t%s" %(chrom[i], start[i], end[i], strand[i], sequence_motif[i])
	print motif_coordinates
	count +=1
print "Total site :", count

######
######
### For motif analysis using fastinterval and bedtools logic, so this is zero-based so the end element is exclusive:
bin_size = 1
upstream = 50
downstream = 50
bins = range(-upstream, (downstream+1), 1)


test_data = [] 
analysis_dict = {}

#for i in range(5):
#    print analysis_dict["chr1"][i][3]

final_data = list()
final_header = ["chrom", "start", "end", "bin_start", "mid_point", "bin_end", "sequence", "motif_start", "motif_end"]

for i in range(len(final_header)):
	final_data.append(list())

coordinate_file = "/home/surya/surya/for_gata_parse/gata_motif.bed"
with open(coordinate_file, "w") as outfile:
	print "chrom\tmid_point\tstart\tend\tbin_start\tbin_end\tsequence\tmotif_start\tmotif_end"
	for i in range(len(chrom)):
		chrome = chrom[i]
		print chrome
		motif_start = int(start[i])
		motif_end = int(end[i])
		mid_point = (motif_start + motif_end) / 2
		sequence = sequence_motif[i]

		for each in bins :
			bin_start = each
			bin_end = (each + bin_size)
			#For fastinterval and bedtools case:
			#bin_end = (each + bin_size)
			chrom_start = mid_point + bin_start
			chrom_end = mid_point + bin_end

			line_coords = [chrome, str(chrom_start), str(chrom_end), str(bin_start), str(mid_point), str(bin_end), sequence, str(motif_start), str(motif_end)]
			test_data.append("\t".join(line_coords))

			if chrome in analysis_dict:
				analysis_dict[chrome].append((chrome, int(chrom_start), int(chrom_end), "\t".join(line_coords)))
			else:
				analysis_dict[chrome] = [(chrome, int(chrom_start), int(chrom_end), "\t".join(line_coords))]

				
			for i,item in enumerate(line_coords):
				final_data[i].append(item)

			print_coords = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(chrome, chrom_start, chrom_end, bin_start, mid_point, bin_end, sequence, motif_start, motif_end) 
			print print_coords
			outfile.write(print_coords + "\n")

	final_zip = zip(final_header, final_data)
	final_dict = dict(final_zip)


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
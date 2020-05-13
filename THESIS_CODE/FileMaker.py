'''
This program just takes the updated RNA seq text file with the uniprot data manually added
and makes an updated gene_locus -> common name txt file that contains the STRING, BioCyc and Uniprot
data.

10/17/19
'''
def main():
	file_name = "UNIPROT_RNA_ALL.txt"
	Output_File = "NameDict_w_Uniprot.txt"
	locus_name_Dict, locus_list = File_Reader(file_name)
	File_Writer(locus_list, locus_name_Dict, Output_File)
	return

def File_Reader(file):
	f = open(file, 'r')
	print("opening ", file)
	locus_name_Dict = {}
	locus_list = []
	for l in f:
		line = str(l) #line is equal to string
		itemsinline = line.split()
		locus = itemsinline[0]
		name = itemsinline[1]
		locus_list.append(locus)
		locus_name_Dict[locus] = name
	return locus_name_Dict, locus_list

def File_Writer(locus_list, locus_name_Dict, Output_File):
	f = open(Output_File,'w')
	print("opening ", Output_File)
	for locus in locus_list:
		name = locus_name_Dict[locus]
		f.write(locus + '\t' + name + '\n')
	print("Done Printing to File")
	return

main()
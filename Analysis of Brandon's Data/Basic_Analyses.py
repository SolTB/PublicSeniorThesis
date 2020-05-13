'''
This is just a really basic file to figure out the top up/down-regulated genes in the original file
'''

def main():
	file = "UNIPROT_RNA_ALL.txt"
	UpFile = "Top_Up_Genes.txt"
	DownFile = "Top_Down_Genes.txt"

	locus_name_Dict, locus_list, locus_log2_Dict, log2_list = File_Reader(file)
	Analyzer(locus_name_Dict, locus_list, locus_log2_Dict, log2_list, UpFile, DownFile)
	return

#this just reads the RNASeq file and makes a bunch of data structures that will be useful later on.
def File_Reader(file):
	f = open(file, 'r')
	print("opening ", file)

	locus_name_Dict = {} #locus->common name
	locus_log2_Dict = {} #locus->expression
	locus_list = [] #list of all gene locuses
	log2_list = []

	for l in f:
		line = str(l) #line is equal to string
		itemsinline = line.split()
		locus = itemsinline[0]
		name = itemsinline[1]
		log2_str = itemsinline[3]
		if log2_str != 'NA':
			log2 = float(log2_str)
			sublist = [log2,locus, name]
			log2_list.append(sublist) #only added to the sublist if it has a log 2 value
			locus_log2_Dict[locus] = log2
		locus_list.append(locus)
		locus_name_Dict[locus] = name

	return locus_name_Dict, locus_list, locus_log2_Dict, log2_list

#Input: locus->name dictionary, list of gene locuses, locus->log2 dictionary, list of lists containing gene name/locus and log2fold change
# and names of the ouput files
#Output: Up_Genes.txt (file containing gene locus, common name and log2 fold change of top 100 up regulated)
#and Down_Genes.txt (file containing gene locus, common name and log2 fold change of top 100 down regulated)
def Analyzer(locus_name_Dict, locus_list, locus_log2_Dict, log2_list, UpFile, DownFile):
	negative_list = []
	positive_list = []
	for sublist in log2_list:
		if sublist[0] >=0:
			positive_list.append(sublist)
		else:
			negative_list.append(sublist)

	sorted(positive_list, key= lambda x:x[0])
	top100 = positive_list[0:99]

	Up = open(UpFile, 'w')
	print("opening ", UpFile)
	for sublist in top100:
		Up.write(str(sublist[0]) + '\t' + str(sublist[1]) + '\t' + str(sublist[2])  + '\n')
	Up.close()

	sorted(negative_list, key= lambda x:x[0])
	bottom100 = negative_list[0:99]

	Down = open(DownFile, 'w')
	print("opening ", DownFile)
	for sublist in bottom100:
		Down.write(str(sublist[0]) + '\t' + str(sublist[1]) + '\t' + str(sublist[2])  + '\n')
	Down.close()
	return

main()
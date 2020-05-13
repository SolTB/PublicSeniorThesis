''''
This program takes the original RNA-Seq data and makes a file containing only the 911 genes with names (~16% of the total data) and makes a file for input
into fgsea
Date: 02/10/20
'''

def main():
	InputFile = "TabALL.txt"
	OutputFile = "GSEAInput.txt"
	CommonNameDict = filereader(InputFile)
	OutputMaker(OutputFile, CommonNameDict)
	return


#Input: DESeq output file with common names
#Output: Common Names (NOT EC) -> log2FC dictionary
def filereader(file):
	CommonNameDict = {}
	f=open(file,'r')
	print("Opening " + file)
	for l in f:
		line = str(l)
		itemsinline = line.split()
		name = itemsinline[7]
		l2FC = itemsinline[2]
		if name[0:2] != 'EC' and name[0:2] != 'pO':
			Bigname = name.upper()
			CommonNameDict[Bigname] = l2FC
	return CommonNameDict

#writes output to simple txt file
def OutputMaker(Output_file, CommonNameDict):
	f=open(Output_file,'w')
	print("Writing to file " + Output_file)

	for gene in CommonNameDict:
		f.write(gene + '\t' + CommonNameDict[gene] + '\n')

main()
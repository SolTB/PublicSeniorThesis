
def main():
	InputALL = "EditedFinalResults.txt"
	RankedInput = "RankedGSEAInput.txt"
	Output = "p001Ranked.txt"
	AllGeneList, RankedDict = RankRead(RankedInput)
	NamePadjDict = filereader(InputALL)
	OutputMaker(Output, NamePadjDict, AllGeneList, RankedDict)
	return

def RankRead(file):
	RankedDict = {}
	AllGeneList = [] #list of all genes for which the common name is known

	f=open(file,'r')
	print("Opening " + file)

	for l in f:
		line = str(l)
		itemsinline = line.split('\t')
		name = itemsinline[0]
		log2FC = itemsinline[1]
		log2FC = log2FC.strip('\n')
		if name != "Name":
			AllGeneList.append(name)
			RankedDict[name] = log2FC
	return AllGeneList, RankedDict


def filereader(file):
	NamePadjDict = {} #common name -> padj dict
	f=open(file,'r')
	print("Opening " + file)

	for l in f:
		line = str(l)
		itemsinline = line.split('\t')
		gene_name = itemsinline[7]
		padj = itemsinline[6]
		if padj == "NA" or padj == "P_adj":
			padj = "10" #redefine
		padjnew = float(padj)
		NamePadjDict[gene_name] = padjnew

	return NamePadjDict


def OutputMaker(Output_file, NamePadjDict, AllGeneList, RankedDict):
	f=open(Output_file,'w')
	print("Writing to file " + Output_file)

	count = 0
	for gene in AllGeneList:
		pval = NamePadjDict[gene]
		if pval <= 0.001:
			substring = gene + '\t' + RankedDict[gene]
			f.write(substring + '\n')
			count = count +1
	print("#Genes = ", count)

main()

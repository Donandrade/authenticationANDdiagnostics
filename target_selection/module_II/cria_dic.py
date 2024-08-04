
def dictionary_criate(fasta):
	dic = {}
	for line in fasta.readlines():
		if (line[0] == ">"):
			gi_code = line.strip()
			dic[gi_code] = ""
		else:
			dic[gi_code] = dic[gi_code] + line.strip()
	return dic

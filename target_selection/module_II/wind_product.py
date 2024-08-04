
def janela(fasta_dictionary, nome, wind, tj):
	janela1 ={}
	for i, j in fasta_dictionary.iteritems():
		for a in range(len(j)):
			if len(j[a:1+wind+a]) != tj:
				break
			janela1.setdefault(i, []).append(list(j[a:1+wind+a]))
						
	return janela1

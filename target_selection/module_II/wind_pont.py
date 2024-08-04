def window_pontuation(wind, j_out, t_janela):
        k = 0
        z = 0
        b = 0
#	seq1 = list()
#	seq2 = list()
	for v in wind.values()[k]:
		a = 0
	        z = 0
		x = 0
	        soma1 = 0
	        soma2 = 0
	        pont1 = list()
	        pont2 = list()
		seq1 = list()
		seq2 = list()
		for i in range(len(v)):
			seq1.append(v[a])
			seq2.append(wind.values()[1][b][z])
			if v[a] == "-" or wind.values()[1][b][z] == "-":
				if v[a] == "-":
					pont2.append(3)
					pont1.append(0)		
				elif wind.values()[1][b][z] == "-":
					pont1.append(3)
					pont2.append(0)
			elif v[a]  ==  wind.values()[1][b][z]:
				pont2.append(0)
				pont1.append(0)
	                elif (v[a] == "G" or v[a] == "A" or v[a] == "P" or v[a] == "V" or v[a] =="L" or v[a] == "I" or v[a] == "M") and (wind.values()[1][b][z] == "G" or wind.values()[1][b][z] == "A" or wind.values()[1][b][z] == "P" or wind.values()[1][b][z] == "V" or wind.values()[1][b][z] == "L" or wind.values()[1][b][z] == "I" or wind.values()[1][b][z] == "M"): #"nompol_ali"
	                        pont2.append(1)
	                        pont1.append(1)
	                elif (v[a] == "S" or v[a] == "T" or v[a] == "C" or v[a] == "N" or v[a] == "Q") and (wind.values()[1][b][z] == "S" or wind.values()[1][b][z] == "T" or wind.values()[1][b][z] == "C" or wind.values()[1][b][z] == "N" or wind.values()[1][b][z] == "Q"): #"pol_unch"
	                        pont2.append(1)
	                        pont1.append(1)
	                elif (v[a] == "K" or v[a] == "H" or v[a] == "R") and (wind.values()[1][b][z] == "K" or wind.values()[1][b][z] == "H" or wind.values()[1][b][z] == "R"): #"posit"
	                        pont2.append(1)
	                        pont1.append(1)
	                elif (v[a] == "F" or v[a] == "Y" or v[a] == "W") and (wind.values()[1][b][z] == "F" or wind.values()[1][b][z] == "Y" or wind.values()[1][b][z] == "W"): #"arom"
	                        pont2.append(1)
	                        pont1.append(1)
	                elif (v[a] == "D" or v[a] == "E") and (wind.values()[1][b][z] == "D" or wind.values()[1][b][z] == "E"): #negat
	                        pont2.append(1)
	                        pont1.append(1)
	        
			else:
				pont2.append(2)
				pont1.append(2)
			z += 1
			a += 1
		soma1 = sum(pont1)
		soma2 = sum(pont2)
		pont1 = str(pont1)
		pont2 = str(pont2)
		seq1 = str(seq1)
		seq2 = str(seq2)
#		if((soma1 > 15 and soma1 < 40) and (soma2 > 15 and soma2 < 40)):
		if(soma1 == 45 or soma2 == 45):
			f = open(j_out, 'a')
			f.write((((str(wind.keys()[0]+"\t"+(str(wind.keys()[1]+"\t"+seq1.replace(',', '').replace('\' \'', '').replace('[', '').replace(']', '')+"\t"+str(soma1)+"\t")))))))
			f.write(((seq2.replace(',', '').replace('\' \'', '').replace('[', '').replace(']', '')+"\t"+str(soma2)+"\n")))
			f.close()
		b += 1
		k += 1

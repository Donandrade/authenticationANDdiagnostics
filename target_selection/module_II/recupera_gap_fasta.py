#!/usr/bin/python
#-*- coding: utf-8 -*-

import os

path = "/mnt/c/Users/grafi/Documents/clauida_Embrapa/banana_genomicaComparativa/script/result_proteins"
file_name = list(os.listdir("/mnt/c/Users/grafi/Documents/clauida_Embrapa/banana_genomicaComparativa/script/result_proteins"))

seq = list()
seq1 = list()
a = 0
b = 0

for i in range(len(file_name)):
	file = file_name[i]
	tmp = os.path.join(path, file)
	lines = len(open(tmp).readlines())
	x = 0
	z = 0
	for line in open(tmp).readlines():
		index = list(line.split("\t"))
		x = x+1
		z = z+1
		if(int(index[3]) > 40):
			union = list(index[2])
			if(x != lines):
				seq.append(union[1])#aqui eu vou pegar o primeiro aminoacido de cada janela(nao posso pegar o ultimo. a ultima janela )
			else:
				name = index[0]
				f = open(name.lstrip('>'), 'a')
				f.write(index[0]+"\n"+''.join(seq).replace('-', '')+ str(''.join(union)).replace('\'', '').replace('-', '')+"\n")
		if(int(index[5]) > 40):
			union = list(index[4])
			if(z != lines):		
                        	seq1.append(union[1])#aqui eu vou pegar o primeiro aminoacido de cada janela(nao posso pegar o ultimo. a ultima janela )
                        else:
                                name = index[1]
				f = open(name.lstrip('>'), 'a')
				f.write(index[1]+"\n"+''.join(seq1).replace('-', '')+str(''.join(union)).replace('\'', '').replace('-', '')+"\n")
        seq = list()
	seq1 = list()


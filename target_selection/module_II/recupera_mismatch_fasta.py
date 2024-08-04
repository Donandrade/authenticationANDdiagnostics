#!/usr/bin/python
#-*- coding: utf-8 -*-

import os

path = "/home/es98522/doutorado/cap1/cafe/antigo/widnow_pontuation/mismatches"
file_name = list(os.listdir("/home/es98522/doutorado/cap1/cafe/antigo/widnow_pontuation/mismatches"))

seq = list()
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
		if(index[3] > 17):
#			name = index[0]
#			open(name, 'w') as f
#			f.write(index[0]+"\n"+index[2]+"\n")
			union = list(index[2])
			if(x != lines):
				seq.append(union[1])#aqui eu vou pegar o primeiro aminoacido de cada janela(nao posso pegar o ultimo. a ultima janela )
			else:
				name = index[0]
				f = open(name.lstrip('>'), 'a')
				f.write(index[0]+"\n"+''.join(seq)+ str(''.join(union)).replace('\'', '').replace('-', '')+"\n")
#				print (index[1]+"\n"+''.join(seq)+ str(''.join(union)).replace('\'', '').replace('-', ''))
		elif(index[5] > 17):
			union = list(index[4])
			if(z != lines):		
                        	seq1.append(union[1])#aqui eu vou pegar o primeiro aminoacido de cada janela(nao posso pegar o ultimo. a ultima janela )
                        else:
                                name = index[1]
                                f = open(name.lstrip('>'), 'a')
                                f.write(index[1]+"\n"+''.join(seq1)+str(''.join(union)).replace('\'', '').replace('-', '')+"\n")
			
#				print (index[1]+"\n"+''.join(seq1)+str(''.join(union)).replace('\'', '').replace('-', ''))
        seq = list()
	seq1 = list()


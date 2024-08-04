#!/usr/bin/env python
#-*- coding: utf-8 -*-

#import pandas as pd
import wind_product as wp
import cria_dic as cd
import contaGapAlinhamento as alinha
import subprocess
from subprocess import call
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import re

if (len(sys.argv) <= 4):
	print("Usage: python main.py goodProteins.fasta ortolog_1x1.out alinhamento.out <window_size>")
	exit()

#fasta = open(sys.argv[1])
orthologue = open(sys.argv[2])

dic = {}
identifier = list()

for gene in SeqIO.parse(sys.argv[1], "fasta"):
    dic[gene.id] = str(gene.seq)



for line1 in orthologue.readlines():

	line1.strip()

	ortho = re.split('\t|\s+', line1)

	name = ortho[0]+"_"+ortho[1]

	name1 = ortho[0]+"_"+ortho[1]

	name2 = ortho[1]+"_"+ortho[0]

	with  open('par.aa', 'w') as f:

		f.write(">"+ortho[0]+"\n"+dic[ortho[0]]+"\n"+">"+ortho[1]+"\n"+dic[ortho[1]]+"\n")

	subprocess.call("/home/emasilva/anaconda3/bin/muscle -in par.aa -out "+name1, shell=True)
    
	alinhamento = name1

	print(alinhamento)

	alinha.conta_gap(alinhamento, sys.argv[4])
 

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

if (len(sys.argv) <= 2):
	print("Usage: python recuperaFasta.py fasta listsID")
	exit()

identificadores = open(sys.argv[2], 'r').readlines()

dic = {}

for gene in SeqIO.parse(sys.argv[1], "fasta"):
    dic[gene.id] = str(gene.seq)

for k, v in dic.items():
    for i in identificadores:
        i = i.strip()
        #print(k, i)
        if(k == i):
#            print(k,v, i)
            print(">"+k+"\n"+v)
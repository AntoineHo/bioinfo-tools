#!/usr/bin/python
# -*- coding: utf-8 -*-

contigs = []

f = open("/media/antoine/DATA3/ricciae/miniasm/toplot.bed", "r")
for line in f :
	s = line.split('\t')
	name = s[0]
	contigs.append(s[0])
f.close()

for contig in contigs :
	covPE = []
	f = open("/media/antoine/DATA3/ricciae/miniasm/miniasmVpe.sorted.bam.cov", "r")
	
	a = f.readline().split('\t')
	while a[0] != contig :
		a = f.readline().split('\t')
		continue

	covPE.append(int(a[2]))

	while 1 :
		s = f.readline().split('\t')
		if s[0] == contig:
			covPE.append(int(s[2]))
		else :
			break
	f.close()

	x = [i for i in range(0, len(covPE))]

	import matplotlib as mpl
	import matplotlib.pyplot as plt
	import numpy as np

	plt.style.use("seaborn")
	fig, ax = plt.subplots(figsize=(20,5))

	ax.plot(x,covPE, label="PE", linewidth=0.1)
	ax.set_ylim(-40, 1500)
	ax.legend()
	ax.set_title(contig)
	plt.savefig(contig + ".png", format="png", dpi=600)

#!/usr/bin/env python

import pandas as pd
import sys
import numpy as np

d=pd.read_table(sys.argv[1])
d.sort_values(sys.argv[2], inplace=True)

MARGIN=500000

class Peak:
	chrom=0
	start=0
	end=0
	def __init__(self, chr, ps):
		self.chrom=chr
		self.start=ps
		self.end=ps

	def add_snp(self, chr, ps):
		if ps < self.start:
			self.start=ps
		if ps > self.end:
			self.end=ps


class peakCollection:
	peaks=[]
	def __init__(self):
		self.peaks=[]

	def check_and_add(self, chr, ps):
		if len(self.peaks)==0:
#			print("Initialised array of peaks with ",chr,":",ps)
			x=Peak(chr, ps)
			self.peaks=[x]
			return
		found=0
		for i, peak in enumerate(self.peaks):
			if (chr == peak.chrom) and (ps > (peak.start - MARGIN)) and (ps < (peak.end + MARGIN)):
#				print("New snp ",chr, ":", ps, " added to peak ", peak.chrom, ":", peak.start)
				self.peaks[i].add_snp(chr, ps)
				found=1
		if not found:
#				print("Added new peak: ",chr, ":", ps)
				self.peaks.append(Peak(chr, ps))
#		print("At end:")
#		self.print()
#		print("\n\n")
		
	def print(self):
		for peak in self.peaks:
			print(str(peak.chrom)+"\t"+str(peak.start)+"\t"+str(peak.end));
	
	def extend(self, TOTAL_LENGTH):
		## extend the region around its center so that it spans TOTAL_LENGTH
		## If it is already TOTAL_LENGTH or larger, add 1/4 TOTAL_LENGTH either side
		for i,peak in enumerate(self.peaks):
			if peak.end - peak.start < TOTAL_LENGTH:
				center=peak.start+ (peak.end - peak.start)/2
				peak.start=int(center - TOTAL_LENGTH/2)
				peak.end=int(center + TOTAL_LENGTH/2)
			else:
				peak.start=int(peak.start-TOTAL_LENGTH/4)
				peak.end=int(peak.end+TOTAL_LENGTH/4)
			self.peaks[i]=peak
	


p=peakCollection()
for index, row in d.iterrows():
#	print(index, row)
	p.check_and_add(row[sys.argv[3]], row[sys.argv[4]])

#p.print()
p.extend(1000000)
#print("lol")
p.print()

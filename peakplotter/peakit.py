#!/usr/bin/env python3
import sys

import pandas as pd


MARGIN = 500000

class Peak:
	chrom = 0
	start = 0
	end = 0
	def __init__(self, chr, ps):
		self.chrom = chr
		self.start = ps
		self.end = ps

	def add_snp(self, chr, ps):
		if ps < self.start:
			self.start = ps
		if ps > self.end:
			self.end = ps


class peakCollection:
	def __init__(self):
		self.peaks = []

	def check_and_add(self, chr, ps):
		if len(self.peaks)==0:
			x = Peak(chr, ps)
			self.peaks = [x]
			return
		
		found = 0
		for i, peak in enumerate(self.peaks):
			if (chr == peak.chrom) and ((peak.start - MARGIN) < ps < (peak.end + MARGIN)):
				self.peaks[i].add_snp(chr, ps)
				found=1
		
		if not found:
			self.peaks.append(Peak(chr, ps))

		
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
			
			if peak.start<0:
				peak.start=0
			self.peaks[i]=peak

def peakit(signals, pvalcol, chrcol, pscol):
	d = pd.read_table(signals)
	d.sort_values(pvalcol, inplace=True)

	p = peakCollection()
	for index, row in d.iterrows():
		p.check_and_add(row[chrcol], row[pscol])

	p.extend(1000000)
	p.print()

def _peakit(signals: pd.DataFrame, pvalcol, chrcol, pscol):
	signals.sort_values(pvalcol, inplace=True)

	p = peakCollection()
	for index, row in d.iterrows():
		p.check_and_add(row[chrcol], row[pscol])

	p.extend(1000000)
	p.print()

if __name__ == '__main__':
	peakit(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

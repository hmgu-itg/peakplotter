#!/usr/bin/env python3
import sys
from operator import attrgetter

import pandas as pd


class Peak:
    def __init__(self, chrom: int, start: int, end: int):
        self.chrom = chrom
        self.start = start
        self.end = end

    def add_snp(self, chrom: int, pos: int):
        """
        Extend the peak region up to the input SNP's position. 
        
        Raises
        ------
        ValueError
            If the input chromosome is different from the peak's chromosome
        """
        if self.chrom != chrom:
            raise ValueError("Input SNP's chromosome different from the peak's chromosome")
        
        if pos < self.start:
            self.start = pos
        if pos > self.end:
            self.end = pos
    
    def __repr__(self):
        return f'Peak({self.chrom}, {self.start}, {self.end})'
    

class PeakCollection(list):
    def __init__(self, MARGIN = 500_000):
        self.MARGIN = MARGIN
        
    def check_and_add(self, chrom, pos):
        if len(self)==0:
            peak = Peak(chrom, pos, pos)
            self.append(peak)
            return

        found = False
        for i, peak in enumerate(self):
            if (chrom == peak.chrom) and ((peak.start - self.MARGIN) < pos < (peak.end + self.MARGIN)):
                self[i].add_snp(chrom, pos)
                found = True

        if not found:
            self.append(Peak(chrom, pos, pos))
            
    def extend_peaks(self, length: int):
        """
        Extend the region around its center so that it spans `length`.
        If it is already `length` or larger, add 1/4 `length` either side
        """
        for i, peak in enumerate(self):
            if peak.end - peak.start < length:
                center = peak.start + (peak.end - peak.start)/2
                peak.start = int(center - length/2)
                peak.end = int(center + length/2)
            else:
                peak.start = int(peak.start - length/4)
                peak.end = int(peak.end + length/4)

            if peak.start < 0:
                peak.start = 0
                
    @property
    def data(self):
        d = list()
        for peak in self:
            assert peak.start < peak.end
            d.append([peak.chrom, peak.start, peak.end])
        d = pd.DataFrame(d, columns = ['chrom', 'start', 'end']).sort_values(['chrom', 'start']).reset_index(drop=True)
        return d

    def merge(self):
        self.sort()
        merged_collection = PeakCollection()
        curr_peak = None
        for i in range(len(self)):
            if curr_peak is None:
                curr_peak = self[i]
                continue
            
            next_peak = self[i]
            if curr_peak.chrom != next_peak.chrom:
                merged_collection.append(curr_peak)
                curr_peak = next_peak
                continue
            
            curr_peak_range = set(range(curr_peak.start, curr_peak.end + 1))
            next_peak_range = set(range(next_peak.start, next_peak.end + 1))
            
            if curr_peak_range.intersection(next_peak_range):
                curr_peak = Peak(curr_peak.chrom, min(curr_peak.start, next_peak.start), max(curr_peak.end, next_peak.end))
            else:
                merged_collection.append(curr_peak)
                merged_collection.append(next_peak)
        else:
            if curr_peak is not None:
                merged_collection.append(curr_peak)
        
        # Replace 
        self[:] = merged_collection

    def __eq__(self, other):
        if not isinstance(other, PeakCollection):
            return False
        if len(self) != len(other):
            return False
        if str(self) != str(other):
            return False
        return True

    def sort(self):
        self[:] = sorted(self, key = attrgetter('chrom', 'start', 'end'))

        
def peakit(signals: pd.DataFrame, pval_col: str, chr_col: str, pos_col: str) -> PeakCollection:
    sorted_signals = signals.sort_values(pval_col)
    peaks = PeakCollection()
    for index, row in sorted_signals.iterrows():
        peaks.check_and_add(row[chr_col], row[pos_col])
    peaks.extend_peaks(1_000_000)
    return peaks


if __name__ == '__main__':
    peakit(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

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
        if self.data.empty:
            return 

        input_list = list()
        for chrom in set(self.data['chrom']):
            df = self.data.loc[self.data['chrom']==chrom].reset_index(drop=True)
            curr_peak = None
            for i, row in df.iterrows():
                if curr_peak is None:
                    curr_peak = Peak(row['chrom'], row['start'], row['end'])
                    continue
                x = range(curr_peak.start, curr_peak.end)
                y = range(row['start'], row['end'])
                overlap = bool(range(max(x[0], y[0]), min(x[-1], y[-1])+1))
                if overlap:
                    curr_peak = Peak(chrom, min(x[0], y[0]), max(x[-1], y[-1])+1)
                else: 
                    input_list.append(curr_peak)
                    curr_peak = None
                    # We don't want to accidentally discard the last row
                    if i == df.shape[0]-1:
                        curr_peak = Peak(row['chrom'], row['start'], row['end'])
            else:
                input_list.append(curr_peak)
            
        # Replace 
        self[:] = input_list
        self.sort()

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

    @classmethod
    def from_list(cls, input: list, MARGIN = 500_000):
        if not all([isinstance(i, Peak) for i in input]):
            raise ValueError("Input list has object other than 'Peak'.")
        instance = cls(MARGIN)
        instance[:] = input
        
        return instance

        
def peakit(signals: pd.DataFrame, pval_col: str, chr_col: str, pos_col: str, flank: int) -> PeakCollection:
    sorted_signals = signals.sort_values(pval_col)
    peaks = PeakCollection()
    for index, row in sorted_signals.iterrows():
        peaks.check_and_add(row[chr_col], row[pos_col])
    peaks.extend_peaks(flank * 2)
    return peaks


if __name__ == '__main__':
    peakit(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

from peakplotter.peakit import Peak, PeakCollection

def test_Peak_add_snp():
    peak = Peak(1, 500_000, 500_000)
    peak.add_snp(1, 499_999)
    peak.add_snp(1, 500_001)
    assert peak.chrom == 1
    assert peak.start == 499_999
    assert peak.end == 500_001
    
def test_PeakCollection():
    peaks = PeakCollection()
    
    peaks.check_and_add(1, 100)
    peaks.check_and_add(2, 100)
    
    peaks.check_and_add(1, 99)
    peaks.check_and_add(1, 101)
    
    peaks.check_and_add(2, 99)
    peaks.check_and_add(2, 101)
    
    peak1 = peaks[0]
    assert peak1.chrom == 1
    assert peak1.start == 99
    assert peak1.end == 101
    
    peak2 = peaks[1]
    assert peak2.chrom == 2
    assert peak2.start == 99
    assert peak2.end == 101
    
    assert len(peaks) == 2

def test_PeakCollection_merge():
    peaks = PeakCollection()
    peaks[:] = [
        Peak(1, 100, 200),
        Peak(2, 100, 200),
        Peak(1, 150, 300),
    ]

    expected = PeakCollection()
    expected[:] = [
        Peak(1, 100, 300),
        Peak(2, 100, 200),
    ]

    peaks.merge()

    assert peaks == expected
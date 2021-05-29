from pathlib import Path
from peakplotter.tools import Plink

base = Path(__file__).parent

exclude_dir = base.joinpath('exclude')
merge_dir = base.joinpath('merge')


def setup():
    return Plink(1000)

# TODO: Add cleanup
def test_exclude_all_variants_raise_non_zero():
    plink = setup()
    bfile = exclude_dir.joinpath('test')
    excludelist = exclude_dir.joinpath('excludelist')
    out = exclude_dir.joinpath('test_output')

    ps = plink.exclude(bfile, excludelist, out)

    assert ps.returncode!=0, 'Excluding all variants returned exit code 0'
    
from pathlib import Path
from peakplotter.tools import Plink, plink_exclude_across_bfiles

base = Path(__file__).parent

exclude_dir = base.joinpath('exclude')
merge_dir = base.joinpath('merge')


def setup():
    return Plink(1000)

# TODO: Add cleanup
def test_exclude_all_variants_raise_non_zero():
    plink = setup()
    bfile = exclude_dir.joinpath('cohortA')
    excludelist = exclude_dir.joinpath('excludelist')
    out = exclude_dir.joinpath('exclude_all_test')

    ps = plink.exclude(bfile, excludelist, out)

    assert ps.returncode!=0, 'Excluding all variants returned exit code 0'

# TODO: Add cleanup
def test_plink_exclude_across_bfiles():
    """
    cohortA has 2 variants and cohortB has 3 variants.
    Excludelist has 2 variants, which results in the exclusion of
    all variants in cohortA and 2 variants of cohortB.
    The output list of `plink_exclude_across_bfiles` function
    is therefore expected to be a `[Path(f'{cohortB_bfile}.excluded')]`.
    """
    plink = setup()
    cohortA_bfile = exclude_dir.joinpath('cohortA')
    cohortB_bfile = exclude_dir.joinpath('cohortB')
    excludelist = exclude_dir.joinpath('excludelist')

    bfiles = [cohortA_bfile, cohortB_bfile]
    output = plink_exclude_across_bfiles(plink, bfiles, excludelist)
    expected = [Path(f'{cohortB_bfile}.excluded')]

    assert output == expected
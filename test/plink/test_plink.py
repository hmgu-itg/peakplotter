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
    bfile = base.joinpath('cohortA')
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
    cohortA_bfile = base.joinpath('cohortA')
    cohortB_bfile = base.joinpath('cohortB')
    excludelist = exclude_dir.joinpath('excludelist')

    bfiles = [cohortA_bfile, cohortB_bfile]
    output = plink_exclude_across_bfiles(plink, bfiles, excludelist)
    expected = [Path(f'{cohortB_bfile}.excluded')]

    assert output == expected


def test_plink_merge_one_cohort():
    """
    Check that plink just copies the binary file when executing
    merge command with one bfile in mergelist
    """
    plink = setup()
    mergelist = merge_dir.joinpath('mergelist_one_cohort.txt')
    out = merge_dir.joinpath('merge_one_cohort')
    output = plink.merge(mergelist, out)
    
    assert output.returncode == 0, f'Merging one cohort returned exit code {output.returncode}'

def test_plink_merge_two_cohorts():
    """
    Check that plink just copies the binary file when executing
    merge command with one bfile in mergelist
    """
    plink = setup()
    mergelist = merge_dir.joinpath('mergelist_two_cohorts.txt')
    out = merge_dir.joinpath('merge_two_cohorts')
    output = plink.merge(mergelist, out)
    
    assert output.returncode == 0, f'Merging one cohort returned exit code {output.returncode}'


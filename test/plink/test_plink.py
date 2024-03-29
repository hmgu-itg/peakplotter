from pathlib import Path
from peakplotter.tools import Plink, plink_exclude_across_bfiles

base_dir = Path(__file__).parent
# base = Path().absolute()
exclude_dir = base_dir.joinpath('exclude')
merge_dir = base_dir.joinpath('merge')


def setup():
    return Plink(1000)

# TODO: Add cleanup
def test_exclude_all_variants_raise_non_zero():
    plink = setup()
    bfile = base_dir.joinpath('cohortA')
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
    cohortA_bfile = base_dir.joinpath('cohortA')
    cohortB_bfile = base_dir.joinpath('cohortB')
    excludelist = exclude_dir.joinpath('excludelist')

    bfiles = [cohortA_bfile, cohortB_bfile]
    output = plink_exclude_across_bfiles(plink, bfiles, excludelist)
    expected = [Path(f'{cohortB_bfile}.excluded')]

    assert output == expected

def test_plink_exclude_genotypes_from_region_with_no_variants():
    plink = setup()
    cohortA_bfile = base_dir.joinpath('cohortA')
    ps = plink.extract_genotypes(cohortA_bfile, 1, 300, 400, 'test_extract_output')
    
    assert ps.returncode == 12, 'Exit code of extract_genotypes is not 12'

    

def test_plink_merge_one_cohort():
    """
    Check that plink just copies the binary file when executing
    merge command with one bfile in mergelist
    """
    plink = setup()
    mergelist = base_dir.joinpath('mergelist_one_cohort.txt')
    with open(mergelist, 'w') as f:
        f.write(f'{base_dir.joinpath("cohortA")}\n')
    out = base_dir.joinpath('merge_one_cohort')
    output = plink.merge(mergelist, out)
    
    assert output.returncode == 0, f'Merging one cohort returned exit code {output.returncode}'

def test_plink_merge_two_cohorts():
    """
    Check that plink just copies the binary file when executing
    merge command with one bfile in mergelist
    """
    plink = setup()
    mergelist = base_dir.joinpath('mergelist_two_cohorts.txt')
    with open(mergelist, 'w') as f:
        f.write(f'{base_dir.joinpath("cohortA")}\n')
        f.write(f'{base_dir.joinpath("cohortB")}\n')
    out = base_dir.joinpath('merge_two_cohorts')
    output = plink.merge(mergelist, out)
    
    assert output.returncode == 0, f'Merging one cohort returned exit code {output.returncode}'


# PeakPlotter : automatically annotate hits from genome-wide association results

PeakPlotter takes away the annoying task of running regional association plots and annotating variants for your association studies results. It is compatible with sequencing as well as GWAS data. It is compatible with any format (GEMMA, SNPTEST, Bolt-LMM...) that produces the relevant columns: chromosome, position, unique ID, P-value, reference and non-reference alleles.

## Install

After installing the prerequisites (see below), clone the repository and install using `pip`.
```bash
git clone https://github.com/hmgu-itg/peakplotter.git

cd peakplotter

python3 -m pip install .

peakplotter-data-setup # This only needs to be run once

peakplotter --help
# or 
python3 -m peakplotter --help
```

A `Singularity` definition file is also available in the repository if you wish to build a container to use `peakplotter`.

## Prerequisites
PeakPlotter has has non-python dependencies.  
In order to run PeakPlotter you need to install the following tools and add the executables to your `PATH`:
* Plink 1.9 or newer ([available here](https://www.cog-genomics.org/plink2/index))
* LocusZoom Standalone 1.4 or newer ([available here](http://genome.sph.umich.edu/wiki/LocusZoom_Standalone))
* Tabix ([available here](https://github.com/samtools/htslib))
* Moreutils (for `sponge`)

PeakPlotter will throw a `MissingExecutableError` if you have any of the above tools missing in your `PATH` environment variable.  
Add the necessary tools to your `PATH` like below:  
```bash
export PATH=/path/to/locuszoom:/path/to/plink:$PATH
```

If you want to make these changes permanent, do:
```bash
echo 'export PATH=/path/to/locuszoom:/path/to/plink:$PATH' >> ~/.bashrc
```

## Usage
```bash
$ peakplotter --help
Usage: peakplotter [OPTIONS]

  PeakPlotter

Options:
  -a, --assoc-file FILE    Path to the association file. It can be gzipped,
                           provided that it bears the .gz extension. Its first
                           line must be a header, coherent with the name
                           arguments below. It must be tab-separated, bgzipped
                           and tabixed (tabix is available as part of
                           bcftools)  [required]
  -f, --bfiles TEXT        Binary PLINK (.bed/.bim/.fam) file base name. This
                           should contain the genotypes for at least all the
                           variants in the assoc_file, but it can contain
                           more. Please note that this is the base name,
                           without the .bed/.bim/.fam extension.  [required]
  -o, --out DIRECTORY      Output directory to store all output files.
                           [required]
  -chr, --chr-col TEXT     Name of the column for chromosome names.
                           [required]
  -ps, --pos-col TEXT      Name of the column for chromosomal position.
                           [required]
  -rs, --rs-col TEXT       Name of the column for unique SNP ids (RS-id or
                           chr:pos).  [required]
  -p, --pval-col TEXT      Name of the column for p-values.  [required]
  -a1, --a1-col TEXT       Name of the column for reference or major allele
                           (used for predicting consequence).  [required]
  -a2, --a2-col TEXT       Name of the column for alternate or minor allele.
                           [required]
  -maf, --maf-col TEXT     Name of the column for non-reference or minor
                           allele frequency.  [required]
  -b, --build INTEGER      Assembly build (37 or 38)  [default: 38]
  -s, --signif FLOAT       The significance level above which to declare a
                           variant significant. Scientific notation (such as
                           5e-8) is fine.
  -bp, --flank-bp INTEGER  Flanking size in base pairs for drawing plots
                           (defaults to 500kb, i.e. 1Mbp plots) around lead
                           SNPs.
  --overwrite              Overwrite output directory if it already exists.
  --help                   Show this message and exit.
``` 

## Testing
Run `pytest` at the root of the repository to run the testsuite.

There aren't a lot of tests right now, and this is a work in progress. If you encounter any bugs, please raise an issue at the [issue page](https://github.com/hmgu-itg/peakplotter/issues).

# PeakPlotter : automatically annotate hits from genome-wide association results

PeakPlotter takes away the annoying task of running regional association plots and annotating variants for your association studies results. It is compatible with sequencing as well as GWAS data. It is compatible with any format (GEMMA, SNPTEST, Bolt-LMM...) that produces the relevant columns: chromosome, position, unique ID, P-value, reference and non-reference alleles.

## Install
Clone the repository and install using `pip`.
```bash
git clone https://github.com/hmgu-itg/peakplotter.git

cd peakplotter

python3 -m pip install .

peakplotter --help
# or 
python3 -m peakplotter --help
```

## Prerequisites
PeakPlotter also have non-python dependencies.  
In order to run PeakPlotter you need to install the following tools and add the executables to your `PATH`:
* Plink 1.9 or newer ([available here](https://www.cog-genomics.org/plink2/index))
* LocusZoom Standalone 1.3 or newer ([available here](http://genome.sph.umich.edu/wiki/LocusZoom_Standalone))
* BedTools ([available here](http://bedtools.readthedocs.io/en/latest/))
* Tabix ([available here](https://github.com/samtools/htslib))
* Coreutils (for `sponge`)

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
$ python3 -m peakplotter --help
Usage: python -m peakplotter [OPTIONS]

  PeakPlotter

Options:
  -a, --assoc-file FILE    Path to the association file. It can be gzipped,
                           provided that it bears the .gz extension. Its first
                           line must be a header, coherent with the name
                           arguments below. It must be tab-separated, bgzipped
                           and tabixed (tabix is available as part of
                           bcftools)
  -f, --bfiles FILE        Binary PLINK (.bed/.bim/.fam) file base name. This
                           should contain the genotypes for at least all the
                           variants in the assoc_file, but it can contain
                           more. Please note that this is the base name,
                           without the .bed/.bim/.fam extension.
  -s, --signif FLOAT       The significance level above which to declare a
                           variant significant. Scientific notation (such as
                           5e-8) is fine.
  -chr, --chr-col TEXT     Name of the column for chromosome names.
  -ps, --pos-col TEXT      Name of the column for chromosomal position.
  -rs, --rs-col TEXT       Name of the column for unique SNP ids (RS-id or
                           chr:pos).
  -p, --pval-col TEXT      Name of the column for p-values.
  -a1, --a1-col TEXT       Name of the column for reference or major allele
                           (used for predicting consequence).
  -a2, --a2-col TEXT       Name of the column for alternate or minor allele.
  -maf, --maf-col TEXT     Name of the column for non-reference or minor
                           allele frequency.
  -bp, --flank-bp INTEGER  Flanking size in base pairs for drawing plots
                           (defaults to 500kb, i.e. 1Mbp plots) around lead
                           SNPs.
  --ref-flat FILE          Path to Locuszoom\'s refFlat file. By default,
                           peakplotter finds it for you in the locuszoom
                           files.
  --recomb FILE            Path to Locuszoom\'s recomb_rate file. By default,
                           peakplotter finds it for you in the locuszoom
                           files.
  --overwrite              Overwrite output directory if it already exists.
  --help                   Show this message and exit.
``` 

## Output
Curently, the script outputs results in the current directory. More precisely, it appends extensions to the association results filename, so it is generally good to assume that PeakPlotter is not very good at handling relative, or even absolute paths. It is therefore safest to work in the same directory as the association results file. Paths to the bed file are insensitive to this issue and can be located wherever you want.

## Troubleshooting
At the moment PeakPlotter does not handle errors very well. In particular, it doesn't catch errors thrown by Plink and doesn't stop if something goes wrong at some point in the pipeline. This is work in progress, if you want to report a bug, keep your logs and email [the author](mailto:ag15@sanger.ac.uk).

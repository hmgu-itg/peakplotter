# PeakPlotter : automatically annotate hits from genome-wide association results

PeakPlotter takes away the annoying task of running regional association plots and annotating variants for your association studies results. It is compatible with sequencing as well as GWAS data. It is compatible with any format (GEMMA, SNPTEST, Bolt-LMM...) that produces the relevant columns: chromosome, position, unique ID, P-value, reference and non-reference alleles.

## Prerequisites
In order to run PeakPlotter you need to have working copies of the following tools installed:
* Plink 1.9 or newer ([available here](https://www.cog-genomics.org/plink2/index))
* LocusZoom Standalone 1.3 or newer ([available here](http://genome.sph.umich.edu/wiki/LocusZoom_Standalone))
* BedTools ([available here](http://bedtools.readthedocs.io/en/latest/))

### Sanger users
In order to run PeakPlotter on the Sanger servers, please paste the following in your terminal:

```bash
. /software/hgi/etc/profile.hgi
module add $(module avail 2>&1 | grep '/plink/' | grep latest | sed 's/.latest.//')
module add $(module avail 2>&1 | grep '/bedtools/' | grep latest | sed 's/.latest.//')
module add $(module avail 2>&1 | grep '/perl/' | grep latest | sed 's/.latest.//')
export PERL5LIB=/nfs/users/nfs_a/ag15/perlmod/lib/site_perl/5.20.1:$PERL5LIB
export PATH=/nfs/team144/software/locuszoom-1.3/locuszoom/bin:$PATH
source /nfs/team144/jupyter_p27/bin/activate
```

## Installation
Hopefully, PeakPlotter should work quasi out-of-the-box. It **needs LocusZoom, BedTools and Plink to be in your path**, since it calls `plink`, `bedtools` and `locuszoom` directly from inside the script. This is done like:

```bash
export PATH=/path/to/locuszoom:/path/to/plink:$PATH
```

If you want to make these changes permanent, do:
```bash
echo 'export PATH=/path/to/locuszoom:/path/to/plink:$PATH' >> ~/.bashrc
```

Additionally, you will need to modify two paths in the executable. In your favourite text editor, open `plotpeaks.sh` and change the following lines:
```bash
## LOCUSZOOM DATA PATHS
## (N.B. LZ path has to be in PATH)
REFFLAT="/nfs/team144/software/locuszoom-1.2/data/database/refFlat.txt"
RECOMB="/nfs/team144/software/locuszoom-1.2/data/database/recomb-rate.txt"
```

to reflect your LocusZoom install path. For example, if `/path/to/locuszoom` is your locuszoom install path, then the path to the binaries will be `/path/to/locuszoom/bin`, and the `REFFLAT` and `RECOMB` above should be set to `/path/to/locuszoom/data/database/refFlat.txt` and `/path/to/locuszoom/data/database/recomb-rate.txt`.

### Python libraries
You will need the following Python libraries in order to run the program: `pandas, numpy, bokeh, urllib2, json, requests, asr`.
Most of these packages are available in large Python bundles such as Anaconda. For those that are missing, or if you prefer to use your own Python environment, you can test the presence of a module and install it if it is missing:
```bash
python -c "import bokeh" || pip install bokeh
```

## Syntax

```bash
./plotpeaks.sh [signif] [assoc_file] [chrcol_name] [poscol_name] [idcol_name] [pvalcol_name] [allele1col_name] [allele0col_name] [afcol_name] [bed_file]
```
* **signif** is the significance level above which to declare a variant significant. Scientific notation (such as `5e-8`) is fine.
* **assoc_file** is the association file. It can be gzipped, provided that it bears the `.gz` extension. Its first line must be a header, coherent with the name arguments below. It must be **tab-separated**.
* **chrcol_name** : name of the column for chromosome names.
* **poscol_name** : name of the column for chromosomal position.
* **idcol_name** : name of the column for unique SNP ids (RS-id or chr:pos).
* **pvalcol_name** : name of the column for p-values.
* **allele1col_name** : name of the column for reference or major allele (used for predicting consequence).
* **allele0col_name** : name of the column for alternate or minor allele.
* **afcol_name** : name of the column for non-reference or minor allele frequency.
* **bed_file** : BED file base name. This should contain the genotypes for at least all the variants in the **assoc_file**, but it can contain more. Please note that this is the base name, without the `.bed/.bim/.fam` extension.
* **region_flank (optional)** flanking size in base pairs for drawing plots (defaults to 500kb, i.e. 1Mbp plots) around lead SNPs.

## Genome build
At the moment, PeakPlotter is compatible with GRCh37(hg19) and GRCh38. It defaults to the latter, and can be forced to run on b37 by **prefixing all arguments** by `b37`, as in :
```bash
./plotpeaks.sh b37 [signif] [assoc_file] [...]
``` 

## Output

Curently, the script outputs results in the current directory. More precisely, it appends extensions to the association results filename, so it is generally good to assume that PeakPlotter is not very good at handling relative, or even absolute paths. It is therefore safest to work in the same directory as the association results file. Paths to the bed file are insensitive to this issue and can be located wherever you want.

## Troubleshooting
At the moment PeakPlotter does not handle errors very well. In particular, it doesn't catch errors thrown by Plink and doesn't stop if something goes wrong at some point in the pipeline. This is work in progress, if you want to report a bug, keep your logs and email [the author](mailto:ag15@sanger.ac.uk).

# PeakPlotter : automatically annotate hits from genome-wide association results

PeakPlotter takes away the annoying task of running regional association plots and annotating variants for your association studies results. It is compatible with sequencing as well as GWAS data. It is compatible with any format (GEMMA, SNPTEST, Bolt-LMM...) that produces the relevant columns: chromosome, position, unique ID, P-value, reference and non-reference alleles.

## Prerequisites
In order to run PeakPlotter you need to have working copies of the following tools installed:
* Plink 1.9 or newer ([available here](https://www.cog-genomics.org/plink2/index))
* LocusZoom Standalone 1.3 or newer ([available here](http://genome.sph.umich.edu/wiki/LocusZoom_Standalone))

## Installation
Hopefully, PeakPlotter should work quasi out-of-the-box. It **needs LocusZoom and Plink to be in your path**, since it calls `plink` and `locuszoom` directly from inside the script. This is done like:

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

## Syntax

```bash
./plotpeaks.sh [signif] [assoc_file] [chrcol_name] [poscol_name] [idcol_name] [pvalcol_name] [allele1col_name] [allele0col_name] [afcol_name] [bed_file]
```
* **signif** is the significance level above which to declare a variant significant. Scientific notation (such as `5e-8`) is fine.
* **assoc_file** is the association file. It can be gzipped, provided that it bears the `.gz` extension. Its first line must be a header, coherent with the name arguments below.
* **chrcol_name** : name of the column for chromosome names.
* **poscol_name** : name of the column for chromosomal position.
* **idcol_name** : name of the column for chromosomal position.
* **pvalcol_name** : name of the column for chromosomal position.
* **allele1col_name** : name of the column for chromosomal position.
* **allele0col_name** : name of the column for chromosomal position.
* **afcol_name** : name of the column for chromosomal position.
* **bed_file** : BED file base name. This should contain the genotypes for at least all the variants in the **assoc_file**, but it can contain more. Please note that this is the base name, without the `.bed/.bim/.fam` extension.
* **region_flank (optional)** flanking size in base pairs for drawing plots (defaults to 500kb, i.e. 1Mbp plots) around lead SNPs.

## Genome build
At the moment, PeakPlotter is compatible with GRCh37(hg19) and GRCh38. It defaults to the latter, and can be forced to run on b37 by **prefixing all arguments** by `b37`, as in :
```bash
./plotpeaks.sh b37 [signif] [assoc_file] [...]
``` 

## Troubleshooting
At the moment PeakPlotter does not handle errors very well. In particular, it doesn't catch errors thrown by Plink and doesn't stop if something goes wrong at some point in the pipeline. This is work in progress, if you want to report a bug, keep your logs and email [the author](mailto:ag15@sanger.ac.uk).

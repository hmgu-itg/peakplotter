#!/bin/bash

b37=0

if [[ $1 == "b37" ]]; then
	shift
	echo -e "\n\nBuild 37 requested.\n\n"
	b37=1
else
	echo -e "\n\nBuild 38 requested.\n\n"
	b37=""
fi

signif=$1
assocfile=$2
chrcol=$3
pscol=$4
rscol=$5
pvalcol=$6
files=${10}
a1col=$7
a2col=$8
mafcol=$9
flank_bp=${11}
filelist=$files

## LOCUSZOOM DATA PATHS
## (N.B. LZ path has to be in PATH)
if [ -z "$b37" ]
	then
	REFFLAT="/nfs/team144/software/locuszoom-1.3/locuszoom/data/database/refFlat.b38.txt"
	RECOMB="/nfs/team144/software/locuszoom-1.3/locuszoom/data/database/recomb-rate_GRCh38.txt"
else
        REFFLAT="/nfs/team144/software/locuszoom-1.3/locuszoom/data/database/refFlat.txt"
        RECOMB="/nfs/team144/software/locuszoom-1.3/locuszoom/data/database/recomb-rate.txt"	
fi
## SELF DIR
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo Running from $(pwd), executable in $DIR.

echo
echo
echo $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11}
echo
echo



declare -a files=($(echo $files | tr ',' ' '))
declare -a ids=($(seq 0 $(($(echo ${10} | tr ',' '\n'| wc -l)-1)) | tr '\n' ' '))
if [ -z "$flank_bp" ]; then
	flank_bp=500000
fi

ext_flank_bp=$(( $flank_bp + 100000 ))
flank_kb=$(echo $flank_bp | sed 's/...$//')
ext_flank_kb=$(( $flank_kb + 100 ))
if [ -z $(echo $2 | egrep '.gz$') ] ; then
	cat=cat
else
	cat=zcat
fi

echo $a1col

chrcoli=$(grep -w $chrcol <(paste <(seq 1 $($cat $assocfile | head -n1 | tr '\t' '\n'| wc -l)) <($cat $assocfile | head -n1 | tr '\t' '\n')) | cut -f1)
pscoli=$(grep -w $pscol <(paste <(seq 1 $($cat $assocfile | head -n1 | tr '\t' '\n'| wc -l)) <($cat $assocfile | head -n1 | tr '\t' '\n')) | cut -f1)
rscoli=$(grep -w $rscol <(paste <(seq 1 $($cat $assocfile | head -n1 | tr '\t' '\n'| wc -l)) <($cat $assocfile | head -n1 | tr '\t' '\n')) | cut -f1)
pvalcoli=$(grep -w $pvalcol <(paste <(seq 1 $($cat $assocfile | head -n1 | tr '\t' '\n'| wc -l)) <($cat $assocfile | head -n1 | tr '\t' '\n')) | cut -f1)
a1coli=$(grep -w $a1col <(paste <(seq 1 $($cat $assocfile | head -n1 | tr '\t' '\n'| wc -l)) <($cat $assocfile | head -n1 | tr '\t' '\n')) | cut -f1)
a2coli=$(grep -w $a2col <(paste <(seq 1 $($cat $assocfile | head -n1 | tr '\t' '\n'| wc -l)) <($cat $assocfile | head -n1 | tr '\t' '\n')) | cut -f1)

echo COLUMNS $chrcoli . $pscoli . $rscoli . $pvalcoli . $a1coli . $a2coli

echo "Looking for peaks..."
echo "===================="
echo "(p-value $pvalcol - $pvalcoli; $signif ; $assocfile ; $cat)"
echo
echo
curpeak_p=1
first=1
$cat $assocfile | awk -v p=$pvalcoli -v signif=$signif '$p<signif'| sort -k${chrcoli},${chrcoli}n -k ${pscoli},${pscoli}n  | tr ';' '_' | sed 's/\[b38\]//' > signals
cat <($cat $assocfile | head -1) signals | sponge signals
$DIR/peakit.py signals $pvalcol $chrcol $pscol | sort -k1,1n -k2,2n | bedtools merge -i - | while read chr start end; do 

curpeak_c=$chr
curpeak_ps=$start
# For each signal:

#Extract the assoc file with tabix, assigning chr:pos ids if not already present (or rsids)
tabix $assocfile ${chr}:${start}-$end | awk -v cc=$chrcoli -v pc=$pscoli -v ic=$rscoli '{if(NR>1 && $ic!~/\:/ && $ic!~/rs/){$ic=$cc":"$pc}print}' | sed 's/\[b38\]//;s/\[b37\]//'> peakdata

# Create LocusZoom dB
cat  <(echo -e "snp\tchr\tpos") <(awk -v rs=$rscoli -v chr=$chrcoli -v ps=$pscoli 'BEGIN{OFS="\t"} NR>1 {print $rs, $chr, $ps}' peakdata) | tr ' ' '\t' > peakdata.chrpos
awk '{if($1~/:/ && $1 !~/^chr/){$1="chr"$1}print}' peakdata.chrpos | tr ' ' '\t' |sponge peakdata.chrpos
dbmeister.py --db $curpeak_c.$curpeak_ps.db --snp_pos peakdata.chrpos
dbmeister.py --db $curpeak_c.$curpeak_ps.db --refflat $REFFLAT
dbmeister.py --db $curpeak_c.$curpeak_ps.db --recomb_rate $RECOMB

# Correct start if negative
if [ $start -lt 1 ]
then
	sensible_start=1
	echo -e "\n\n\nWARNING\t Negative start position changed to $sensible_start : $curpeak_c $curpeak_ps (1)\n\n\n"
else
	sensible_start=$start
fi

# For each PLINK file, extract genotype data				
for id in "${ids[@]}"
do 

	plink --allow-no-sex --memory 15000 --bfile ${files[$id]} --chr $chr --from-bp $sensible_start --to-bp $end --out $id --make-bed				
	# Use SNPIDs if rsids not available
	awk 'OFS="\t"{if($2!~/:/ && $2!~/^rs/){$2="chr"$1":"$4}print}' $id.bim | tr ';' '_' | sponge $id.bim
	# add chr if not present in SNPID
	awk 'OFS="\t"{if($2~/:/ && $2!~/^chr/){$2="chr"$2}print}' $id.bim | tr ';' '_' | sponge $id.bim
	sed -i 's/\[b38\]//;s/\[b37\]//' $id.bim
done

# If several input files, attempt to merge
seq 0 $(($(echo $filelist | tr ',' '\n'| wc -l)-1)) > mergelist
numfile=$(cat mergelist | wc -l)

if [ "$numfile" -gt "1" ]; then
	echo -e "\n\nAttempting to merge $numfile files...\n\n"
	plink --allow-no-sex --merge-list mergelist --make-bed --out merged --allow-no-sex
			if [ -f merged-merge.missnp ]
				then
				grep -v -w -f merged-merge.missnp peakdata | sponge peakdata
				for id in "${ids[@]}"
					do 
					plink --allow-no-sex --bfile $id --exclude merged-merge.missnp --make-bed --out $id.tmp
					mv $id.tmp.bed $id.bed
					mv $id.tmp.bim $id.bim
					mv $id.tmp.fam $id.fam
				done
				plink --allow-no-sex --merge-list mergelist --make-bed --out merged --allow-no-sex
			fi
else
	echo -e "\n\nOnly one genotype file provided. Moving on...\n\n"
	mv 0.bim merged.bim
	mv 0.bed merged.bed
	mv 0.fam merged.fam
fi

# Correct IDs in peakdata file
cat <(head -1 signals) peakdata > peakdata.header
# Use SNPIDs if rsids not available
awk 'OFS="\t"{if($2!~/:/ && $2!~/^rs/){$2="chr"$1":"$3}print}' peakdata.header | sponge peakdata.header
# add chr if not present in SNPID
awk 'OFS="\t"{if($2~/:/ && $2!~/^chr/){$2="chr"$2}print}' peakdata.header | sponge peakdata.header


# Extract top SNP
refsnp=$(python -c 'import pandas as pd; d=pd.read_table("peakdata.header");id=d["'$rscol'"][d["'$pvalcol'"].idxmin()];print(id if id.startswith("rs") else "chr"+str(d["'$chrcol'"][d["'$pvalcol'"].idxmin()])+":"+str(int(d["'$pscol'"][d["'$pvalcol'"].idxmin()])))')

echo -e "\n\nIn region $chr $start $end, top SNP is $refsnp\n\n"

# Calculate LD
plink --allow-no-sex --bfile merged --r2 --ld-snp $refsnp  --ld-window-kb $ext_flank_kb $flank_kb_ext --ld-window 999999 --ld-window-r2 0 --out merged
cat <(echo -e "snp1\tsnp2\tdprime\trsquare") <(tail -n +2  merged.ld| awk 'OFS=" "{print $3,$6,$7,$7}') |sed 's/\[b38\]//' |sponge merged.ld

# Running LZ (LZ is build aware)
locuszoom --metal peakdata.header --refsnp "$refsnp" --markercol "$rscol" --pvalcol "$pvalcol" --db $chr.$start.db --prefix $chr.$start.$end.500kb --plotonly showAnnot=T showRefsnpAnnot=T annotPch="21,24,24,25,22,22,8,7" rfrows=20 geneFontSize=.4 --ld merged.ld --start=$sensible_start --end=$end --chr=$chr showRecomb=T

# interactive manh expects comma-separated, with a ld column
join --header -1 $rscoli -2 1 <(cat <(head -n1 peakdata.header) <(tail -n+2 peakdata.header | sort -k$rscoli,$rscoli)) <(cat <(echo $rscol ld) <(cut -f2,3 -d' ' merged.ld |sort -k1,1) ) |tr ' ' ','> $chr.$start.$end.peakdata.ld


# Running interactive manhattan
	if [ -z "$b37" ]
	then
		#echo $DIR/scripts/interactive_manh.py $chr.$start.$end.peakdata.ld "$pvalcol" "$pscol" "$rscol" "$mafcol" "$chrcol" "$a2col" "$a1col" b38
		$DIR/scripts/interactive_manh.py $chr.$start.$end.peakdata.ld "$pvalcol" "$pscol" "$rscol" "$mafcol" "$chrcol" "$a2col" "$a1col" b38
	else
		$DIR/scripts/interactive_manh.py $chr.$start.$end.peakdata.ld "$pvalcol" "$pscol" "$rscol" "$mafcol" "$chrcol" "$a2col" "$a1col" b37
	fi
done

#rm cp* merge* peak* 0.* *.db *signal* *.line  

#!/bin/bash
signif=$1
assocfile=$2
chrcol=$3
pscol=$4
rscol=$5
pvalcol=$6
a1col=$7
a2col=$8
mafcol=$9
files=${10}
flank_bp=${11}
REFFLAT=${12}
RECOMB=${13}
build=${14}
work_dir=${15}
filelist=$files
memory=30000


cd $work_dir

## SELF DIR
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo Running from $(pwd), executable in $DIR.

echo
echo
echo $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15}
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
echo $DIR/peakit.py signals $pvalcol $chrcol $pscol
$DIR/peakit.py signals $pvalcol $chrcol $pscol | sort -k1,1n -k2,2n | bedtools merge -i - > peaked

numpeaks=$(cat peaked | wc -l)

current=1

#cat peaked |  while read chr start end; do

while [ "$numpeaks" -ge "$current" ]
do
  read chr start end <<< $(head -n $current peaked | tail -1)
  echo "Treating peak $chr $start $end (peak $current / $numpeaks )"
  current=$(( $current + 1 ))
  curpeak_c=$chr
  curpeak_ps=$start
  # For each signal:

  # Extract the assoc file with tabix, assigning chr:pos ids if not already present (or rsids)
  # tabix $assocfile ${chr}:${start}-$end | grep -v nan | awk -v cc=$chrcoli -v pc=$pscoli -v ic=$rscoli '{if(NR>1 && $ic!~/\:/ && $ic!~/rs/){$ic=$cc":"$pc}print}' | sed 's/\[b38\]//;s/\[b37\]//'> peakdata
  zcat $assocfile | awk \
    -v curpeak_c="$curpeak_c" \
    -v chrcoli="$chrcoli" \
    -v rscoli="$rscoli" \
    -v pscoli="$pscoli" \
    -v start="$start" \
    -v end="$end" \
    'BEGIN{FS="\t";OFS="\t";}{if ($chrcoli==curpeak_c && start<$pscoli && $pscoli<end){$rscoli=$chrcoli":"$pscoli ; print}}' |
    grep -v nan | sed 's/\[b38\]//;s/\[b37\]//' > peakdata
  # Create LocusZoom dB
  cat  <(echo -e "snp\tchr\tpos") <(awk -v rs=$rscoli -v chr=$chrcoli -v ps=$pscoli 'BEGIN{OFS="\t"} NR>1 {print $rs, $chr, $ps}' peakdata) | tr ' ' '\t' > peakdata.chrpos
  awk '{if($1~/:/ && $1 !~/^chr/){$1="chr"$1}print}' peakdata.chrpos | tr ' ' '\t' | sponge peakdata.chrpos
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
    plink --allow-no-sex --memory $memory --bfile ${files[$id]} --chr $chr --from-bp $sensible_start --to-bp $end --out $id --make-bed
    # Use SNPIDs if rsids not available
    awk 'OFS="\t"{if($2!~/:/ && $2!~/^rs/){$2="chr"$1":"$4}print}' $id.bim | tr ';' '_' | sponge $id.bim
    # add chr if not present in SNPID
    awk 'OFS="\t"{if($2~/:/ && $2!~/^chr/){$2="chr"$2}print}' $id.bim | tr ';' '_' | sponge $id.bim
    sed -i 's/\[b38\]//;s/\[b37\]//' $id.bim
  done

  # If several input files, attempt to merge
  seq 0 $(($(echo $filelist | tr ',' '\n'| wc -l)-1)) > mergelist
  numfile=$(cat mergelist | wc -l)

  if [ "$numfile" -gt "1" ]
  then
    echo -e "\n\nAttempting to merge $numfile files...\n\n"
    plink --allow-no-sex --merge-list mergelist --make-bed --out merged
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
  #<<<<<<< HEAD
  #awk 'OFS="\t"{if($'$rscoli'!~/:/ && $'$rscoli'!~/^rs/){$'$rscoli'="chr"$'$chrcoli'":"$'$pscoli'}print}' peakdata.header | sponge peakdata.header
  #=======
  awk 'OFS="\t"{if(NR>1){if($'$rscoli'!~/:/ && $'$rscoli'!~/^rs/){$'$rscoli'="chr"$'$chrcoli'":"$'$rscoli'}}print}' peakdata.header | sponge peakdata.header
  #>>>>>>> 9db6e2348fe89b634ea16b07d4cd8b7c22132c4d
  # add chr if not present in SNPID
  awk 'OFS="\t"{if(NR>1){if($'$rscoli'~/:/ && $'$rscoli'!~/^chr/){$'$rscoli'="chr"$'$rscoli'}}print}' peakdata.header | sponge peakdata.header

  # Extract top SNP
<<<<<<< HEAD
  refsnp=$(python3 -c 'import pandas as pd; d=pd.read_csv("peakdata.header", sep = "\t");id=d["'$rscol'"][d["'$pvalcol'"].idxmin()];print(id if id.startswith("rs") else "chr"+str(d["'$chrcol'"][d["'$pvalcol'"].idxmin()])+":"+str(int(d["'$pscol'"][d["'$pvalcol'"].idxmin()])))')
=======
  refsnp=$(python -c 'import pandas as pd; d=pd.read_table("peakdata.header");id=d["'$rscol'"][d["'$pvalcol'"].idxmin()];print(id if id.startswith("rs") else "chr"+str(d["'$chrcol'"][d["'$pvalcol'"].idxmin()])+":"+str(int(d["'$pscol'"][d["'$pvalcol'"].idxmin()])))')
>>>>>>> master

  echo -e "\n\nIn region $chr $start $end, top SNP is $refsnp\n\n"

  # Calculate LD
  plink --allow-no-sex --bfile merged --r2 --ld-snp $refsnp  --ld-window-kb $ext_flank_kb $flank_kb_ext --ld-window 999999 --ld-window-r2 0 --out merged
  cat <(echo -e "snp1\tsnp2\tdprime\trsquare") <(tail -n +2  merged.ld | awk 'OFS=" "{print $3,$6,$7,$7}') |sed 's/\[b38\]//' | sponge merged.ld
  
  # Running LZ (LZ is build aware)
  locuszoom \
    --build $build \
    --metal peakdata.header \
    --refsnp "$refsnp" \
    --markercol "$rscol" \
    --pvalcol "$pvalcol" \
    --db $chr.$start.db \
    --prefix $chr.$start.$end.500kb \
    --plotonly showAnnot=T showRefsnpAnnot=T annotPch="21,24,24,25,22,22,8,7" rfrows=20 geneFontSize=.4 \
    --ld merged.ld \
    --start=$sensible_start \
    --end=$end \
    --chr=$chr showRecomb=T

  # interactive manh expects comma-separated, with a ld column
  join --header -1 $rscoli -2 1 <(cat <(head -n1 peakdata.header) <(tail -n+2 peakdata.header | sort -k$rscoli,$rscoli)) <(cat <(echo $rscol ld) <(cut -f2,3 -d' ' merged.ld |sort -k1,1) ) |tr ' ' ','> $chr.$start.$end.peakdata.ld

  # Running interactive manhattan

  if [[ $build -eq 37 ]]
  then
    echo $DIR/interactive_manh.py $chr.$start.$end.peakdata.ld "$pvalcol" "$pscol" "$rscol" "$mafcol" "$chrcol" "$a1col" "$a2col" b37
    $DIR/interactive_manh.py $chr.$start.$end.peakdata.ld "$pvalcol" "$pscol" "$rscol" "$mafcol" "$chrcol" "$a1col" "$a2col" b37
  elif [[ $build -eq 38 ]] 
  then
    echo $DIR/interactive_manh.py $chr.$start.$end.peakdata.ld "$pvalcol" "$pscol" "$rscol" "$mafcol" "$chrcol" "$a1col" "$a2col" b38
    $DIR/interactive_manh.py $chr.$start.$end.peakdata.ld "$pvalcol" "$pscol" "$rscol" "$mafcol" "$chrcol" "$a1col" "$a2col" b38
  fi
  echo "Done with peak $chr $start $end."
done

<<<<<<< HEAD
=======
touch done
>>>>>>> master
# rm merge* peak* 0.* *.db *signal* 

if [ -a cp* ] ; then rm cp*; fi
if [ -a *.line ] ; then rm *.line; fi


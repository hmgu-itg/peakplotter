#!/bin/bash

b37=0

if [[ $1 == "b37" ]]; then
	shift
	echo -e "\n\nBuild 37 requested.\n\n"
	b37=1
else
	echo -e "\n\nBuild 38 requested.\n\n"
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
$cat $assocfile | awk -v p=$pvalcoli -v signif=$signif '$p<signif'| sort -k${chrcoli},${chrcoli}n -k ${pscoli},${pscoli}n  | tr ';' '_' | sed 's/\[b38\]//'  | while read line; 
do

echo $line
#done;
#while false; do
# column numbers. These do not need to be updated every time we loop, but it is needed in some cases and time is not precious here.
chrcoli=$(grep -w $chrcol <(paste <(seq 1 $($cat $assocfile | head -n1 | tr '\t' '\n'| wc -l)) <($cat $assocfile | head -n1 | tr '\t' '\n')) | cut -f1)
pscoli=$(grep -w $pscol <(paste <(seq 1 $($cat $assocfile | head -n1 | tr '\t' '\n'| wc -l)) <($cat $assocfile | head -n1 | tr '\t' '\n')) | cut -f1)
rscoli=$(grep -w $rscol <(paste <(seq 1 $($cat $assocfile | head -n1 | tr '\t' '\n'| wc -l)) <($cat $assocfile | head -n1 | tr '\t' '\n')) | cut -f1)
pvalcoli=$(grep -w $pvalcol <(paste <(seq 1 $($cat $assocfile | head -n1 | tr '\t' '\n'| wc -l)) <($cat $assocfile | head -n1 | tr '\t' '\n')) | cut -f1)
a1coli=$(grep -w $a1col <(paste <(seq 1 $($cat $assocfile | head -n1 | tr '\t' '\n'| wc -l)) <($cat $assocfile | head -n1 | tr '\t' '\n')) | cut -f1)
a2coli=$(grep -w $a2col <(paste <(seq 1 $($cat $assocfile | head -n1 | tr '\t' '\n'| wc -l)) <($cat $assocfile | head -n1 | tr '\t' '\n')) | cut -f1)

	# curXXX contains the current peak info
	# curpeak_XXX contains the previous peak info, same as oldXXX
	curps=$(echo $line | awk -v ps=$pscoli '{print $ps}')
	curchr=$(echo $line | awk -v c=$chrcoli '{print $c}')
	curpval=$(echo $line | awk -v ps=$pvalcoli '{print $ps}')
	curs=$(echo $line | awk -v rs=$rscoli '{print $rs}')
	cura1=$(echo $line | awk -v rs=$a1coli '{print $rs}')
	cura2=$(echo $line | awk -v rs=$a2coli '{print $rs}')
	
	curp=1

	if [[ "$curchr" -ne "$oldchr"  ||  "$curps" -gt "$(($oldpos+$flank_bp))" ]]
		then
		if [ $first -ne 1 ]	
			then
			echo -e "Found peak at $oldchr : $oldpos (p=$curpeak_p)"  
			echo
			echo -e "Fetching info from Ensembl..."
			echo "chr"${curpeak_c}:${curpeak_ps}-${curpeak_ps}_${curpeak_a1}_${curpeak_a2}
			echo "chr"${curpeak_c}:${curpeak_ps}-${curpeak_ps}_${curpeak_a1}_${curpeak_a2}
			if [ -z "$b37" ]
			then
				$DIR/scripts/VarAnnot.b38.pl <(echo "chr"${curpeak_c}:${curpeak_ps}-${curpeak_ps}_${curpeak_a1}_${curpeak_a2}) > $curpeak_c.$curpeak_ps.line
			else
				$DIR/scripts/VarAnnot.pl <(echo "chr"${curpeak_c}:${curpeak_ps}-${curpeak_ps}_${curpeak_a1}_${curpeak_a2}) > $curpeak_c.$curpeak_ps.line
			fi

			#perl /nfs/team144/ds26/FunctionalAnnotation/20150505_new_tool/VarAnnot.pl <(echo "chr"${curpeak_c}:${curpeak_ps}-${curpeak_ps}_${curpeak_a1}_${curpeak_a2}) > $curpeak_c.$curpeak_ps.line
			cp gwas_signals.tsv $curpeak_c.$curpeak_ps.signals
			echo -e "\t extracting $flank_bp" region...
			cat <($cat $assocfile | head -n1 | tr '\t' ' ') <(awk -v chr=$chrcoli -v matchr=$curpeak_c -v ps=$pscoli -v matchps=$curpeak_ps -v flank=$flank_bp '$chr==matchr && $ps>(matchps - flank) && $ps<(matchps+flank)' <($cat $assocfile)) | sed 's/ /\t/g;s/\[b38\]//' > peakdata
			awk -v cc=$chrcoli -v pc=$pscoli -v ic=$rscoli '{if(NR>1 && $ic!~/\:/ && $ic!~/rs/){$ic=$cc":"$pc}print}' peakdata > t
			mv t peakdata



			###### UNCOMMENT THIS FOR CLEANUP OF HORIZONTAL LINES ON PLOT DUE TO SNPS WITH IDENTICAL IMPUTATION
			######### If there is no reason to enable this, leave it out as it WILL mess up the data.



			#			/software/bin/Rscript --vanilla ~/association-scripts/cleanup.R
			#			mv peakdata peakdata.unclean
			#			mv peakdata.cleaned peakdata
			cat  <(echo -e "snp\tchr\tpos") <(awk -v rs=$rscoli -v chr=$chrcoli -v ps=$pscoli 'BEGIN{OFS="\t"} NR>1 {print $rs, $chr, $ps}' peakdata) | tr ' ' '\t' > peakdata.chrpos
			dbmeister.py --db $curpeak_c.$curpeak_ps.db --snp_pos peakdata.chrpos
			dbmeister.py --db $curpeak_c.$curpeak_ps.db --refflat $REFFLAT
			dbmeister.py --db $curpeak_c.$curpeak_ps.db --recomb_rate $RECOMB
			for id in "${ids[@]}"
			do 
			        sensible_start=$(($curpeak_ps-$flank_bp))
				if [ $sensible_start -lt 1 ]
				then
				    sensible_start=1
				    echo -e "\n\n\nWARNING\t Negative start position changed to $sensible_start : $curpeak_c $curpeak_ps (1)\n\n\n"
				fi
				plink --memory 15000 --bfile ${files[$id]} --chr $curpeak_c --from-bp $sensible_start --to-bp $(($curpeak_ps+$flank_bp)) --out $id --make-bed				
				awk 'OFS="\t"{if($2!~/:/ && $2!~/rs/){$2="chr"$1":"$4}print}' $id.bim | tr ';' '_'> t
				mv t $id.bim

			done

			seq 0 $(($(echo $filelist | tr ',' '\n'| wc -l)-1)) > mergelist
			plink --merge-list mergelist --make-bed --out merged --allow-no-sex
			if [ -f merged-merge.missnp ]
				then
				grep -v -w -f merged-merge.missnp peakdata > t 
				mv t peakdata
				for id in "${ids[@]}"
					do 
					plink --bfile $id --exclude merged-merge.missnp --make-bed --out $id.tmp
					mv $id.tmp.bed $id.bed
					mv $id.tmp.bim $id.bim
					mv $id.tmp.fam $id.fam
				done
				plink --merge-list mergelist --make-bed --out merged --allow-no-sex
			fi


			awk '{if($1~/\:/ && !($1~/chr/)){$1="chr"$1}print}' peakdata.chrpos | tr '\t' ' '> t
			mv t peakdata.chrpos
			awk -v rspos=$rscoli '{if($rspos~/\:/ && !($rspos~/chr/)){$rspos="chr"$rspos}print}' peakdata | tr '\t' ' '> t
			mv t peakdata
			awk '{if($2~/\:/ && !($2~/chr/)){$2="chr"$2}print}' merged.bim |sed 's/\[b38\]//'> t
			mv t merged.bim
			refsnp=$(echo $curpeak_rs | awk '{if($0~/\:/ && !($0~/chr/)){$0="chr"$0}print}')
			refsnp=$(echo $refsnp | awk -v c=$curpeak_c -v p=$curpeak_ps '{if($0!~/\:/ && $0!~/rs/){$0="chr"c":"p}print}')
			echo 
			echo
			echo $curpeak_rs $refsnp
			echo
			echo
			plink --bfile merged --r2 --ld-snp $refsnp  --ld-window-kb $ext_flank_kb $flank_kb_ext --ld-window 999999 --ld-window-r2 0 --out merged
			cat <(echo -e "snp1\tsnp2\tdprime\trsquare") <(tail -n +2  merged.ld| awk 'OFS=" "{print $3,$6,$7,$7}') |sed 's/\[b38\]//' > t
			cp t $curpeak_c.$curpeak_ps.ld
			cp $curpeak_c.$curpeak_ps.ld merged.ld
			echo -e "\n\nRunning command:=================\n\n"
		        sensible_start=$(($curpeak_ps-$flank_bp))
				if [ $sensible_start -le 0 ]
				then
				    sensible_start=1
				    echo -e "\n\n\nWARNING\t Negative start position changed to $sensible_start  : $curpeak_c $curpeak_ps (2)\n\n\n"
				fi

			echo locuszoom --metal ../peakdata --refsnp "$refsnp" --markercol "$rscol" --pvalcol "$pvalcol" --db ../$curpeak_c.$curpeak_ps.db --prefix $curpeak_c.$curpeak_ps.500kb --plotonly showAnnot=T showRefsnpAnnot=T annotPch="21,24,24,25,22,22,8,7" rfrows=20 geneFontSize=.4 --ld ../$curpeak_c.$curpeak_ps.ld --start=$sensible_start --end=$(($curpeak_ps+$flank_bp)) --chr=$curpeak_c showRecomb=T --delim \' \'
			#mkdir temp && cd temp
			#~/association-scripts-git/wrap.locuszoom.sh ../peakdata $refsnp $rscol $pvalcol ../$curpeak_c.$curpeak_ps.db $curpeak_c.$curpeak_ps.500kb ../$curpeak_c.$curpeak_ps.ld $sensible_start $(($curpeak_ps+$flank_bp)) $curpeak_c 
			# /nfs/team144/software/locuszoom-1.3/locuszoom/bin/locuszoom --metal ../peakdata --refsnp "$refsnp" --markercol "$rscol" --pvalcol "$pvalcol" --db ../$curpeak_c.$curpeak_ps.db --prefix $curpeak_c.$curpeak_ps.500kb --plotonly showAnnot=T showRefsnpAnnot=T annotPch="21,24,24,25,22,22,8,7" rfrows=20 geneFontSize=.4 --ld ../$curpeak_c.$curpeak_ps.ld --start=$sensible_start --end=$(($curpeak_ps+$flank_bp)) --chr=$curpeak_c showRecomb=T --delim ' '
			#cd -
			#mv temp/* . && rm -rf temp
			#echo -e "\n\n\n$line\n\n\n"
			cp peakdata peakdata.$refsnp.bak
			echo $refsnp $curpeak_c $curpeak_ps>> peaks.plotted
			#rm merged.* peakdata* 
			#cat <(echo $(head -n1 peakdata) ld) <(join -1 $rscoli -2 1 <(sort -k$rscoli,$rscoli peakdata) <(cut -f2,3 -d' ' merged.ld |sort -k1,1) | perl -lane 'print join(@F[2..]), "\n";') > peakdata.ld
			echo JOIN1
			join --header -1 $rscoli -2 1 <(cat <(head -n1 peakdata) <(sort -k$rscoli,$rscoli peakdata)) <(cat <(echo $rscol ld) <(cut -f2,3 -d' ' merged.ld |sort -k1,1) ) > peakdata.ld
			#join --header -j 1 <(cat <(head -n1 peakdata) <(sort -k$rscoli,$rscoli peakdata) | perl -lape 'print join($F[%ENV{"rscoli"}], @F), "\n";') <(cat <(echo $rscol ld) <(cut -f2,3 -d' ' merged.ld |sort -k1,1) ) > peakdata.ld
			## Beware rsid is now number one and the columns are messed up.
			mem=$assocfile
			assocfile=peakdata.ld
			chrcoli=$(grep -w $chrcol <(paste <(seq 1 $(cat $assocfile | head -n1 | tr ' ' '\n'| wc -l)) <(cat $assocfile | head -n1 | tr ' ' '\n')) | cut -f1)
			pscoli=$(grep -w $pscol <(paste <(seq 1 $(cat $assocfile | head -n1 | tr ' ' '\n'| wc -l)) <(cat $assocfile | head -n1 | tr ' ' '\n')) | cut -f1)
			rscoli=$(grep -w $rscol <(paste <(seq 1 $(cat $assocfile | head -n1 | tr ' ' '\n'| wc -l)) <(cat $assocfile | head -n1 | tr ' ' '\n')) | cut -f1)
			pvalcoli=$(grep -w $pvalcol <(paste <(seq 1 $(cat $assocfile | head -n1 | tr ' ' '\n'| wc -l)) <(cat $assocfile | head -n1 | tr ' ' '\n')) | cut -f1)
			a1coli=$(grep -w $a1col <(paste <(seq 1 $(cat $assocfile | head -n1 | tr ' ' '\n'| wc -l)) <(cat $assocfile | head -n1 | tr ' ' '\n')) | cut -f1)
			a2coli=$(grep -w $a2col <(paste <(seq 1 $(cat $assocfile | head -n1 | tr ' ' '\n'| wc -l)) <(cat $assocfile | head -n1 | tr ' ' '\n')) | cut -f1)
			#cat <(echo $(head -n1 peakdata.ld) id gene consequence) <(join -1 3 -2 1 -a 1 <(sort -k3,3 peakdata.ld) <(~ag15/scripts/getrsid.pl <(awk '$NF>0.1' peakdata.ld | tail -n+2 | cut -f2-5 -d' ' )| sort -k1,1) | grep -v MARKER)| awk '{if(NF<22){$21="NA";$22="NA";$23="NA"}if(NR>1){un=$1;deux=$2;trois=$3;$1=deux;$2=trois;$3=un;}print}' > peakdata.ld.annotated
			expnumcol=$(($(head -n1 peakdata.ld | tr ' ' '\n' | wc -l)+3))
			echo JOIN2
			if [ -z "$b37" ]
			then
				join -a 1 -1 $pscoli -2 1 --header <(cat <(head -n1 peakdata.ld) <(tail -n+2 peakdata.ld | sort -k$pscoli,$pscoli)) <( cat <(echo $pscol id gene consequence) <($DIR/scripts/getrsid.b38.pl <(awk '$NF>0.1' peakdata.ld | tail -n+2 | awk -v chr=$chrcoli -v ps=$pscoli -v a1=$a1coli -v a2=$a2coli '{print $chr,$ps,$a1,$a2}') | sort -k1,1))  |awk -v en=$expnumcol '{if(NF<en){$(en-2)="NA";$(en-1)="NA";$en="NA"}print}'> peakdata.ld.annotated
				$DIR/scripts/get_phenotype.b38.pl <(cut -f$(($expnumcol-2)) -d' ' peakdata.ld.annotated | grep -v -e NA -e novel -e id) | sed 's/ /&/g'| sed 's/\t/\t"/;s/$/"/' > phenotypes
			else
				join -a 1 -1 $pscoli -2 1 --header <(cat <(head -n1 peakdata.ld) <(tail -n+2 peakdata.ld | sort -k$pscoli,$pscoli)) <( cat <(echo $pscol id gene consequence) <($DIR/scripts/getrsid.b38.pl b37 <(awk '$NF>0.1' peakdata.ld | tail -n+2 | awk -v chr=$chrcoli -v ps=$pscoli -v a1=$a1coli -v a2=$a2coli '{print $chr,$ps,$a1,$a2}') | sort -k1,1))  |awk -v en=$expnumcol '{if(NF<en){$(en-2)="NA";$(en-1)="NA";$en="NA"}print}'> peakdata.ld.annotated
				$DIR/scripts/get_phenotype.b38.pl b37 <(cut -f$(($expnumcol-2)) -d' ' peakdata.ld.annotated | grep -v -e NA -e novel -e id) | sed 's/ /&/g'| sed 's/\t/\t"/;s/$/"/' > phenotypes
			fi

			
			
			#cat <(echo $(head -n1 peakdata.ld.annotated) assoc) <(join -1 21 -2 1 -a 1 <(tail -n+2 peakdata.ld.annotated | sort -k21,21) <(sort -k1,1 phenotypes) |awk '{if(NF!=24){$24="NA"}print}'| perl -lane 'print join(" ", @F[1..20], $F[0],@F[21..24]);') |sed 's/ /,/g;s/_/ /g'> $curpeak_c.$curpeak_ps.peakdata.ld.annotated.assoc
			echo JOIN3
			join --header -a 1 -1 $(($expnumcol-2)) -2 1 <(cat <(head -n1 peakdata.ld.annotated) <(tail -n+2 peakdata.ld.annotated|  sort -k$(($expnumcol-2)),$(($expnumcol-2)))) <( cat <(echo id assoc) <(sort -k1,1 phenotypes)) |awk -v field=$((expnumcol+2)) '{if(NF!=field){$field="NA"}print}' |sed 's/ /,/g;s/&/ /g'> $curpeak_c.$curpeak_ps.peakdata.ld.annotated.assoc

			if [ -z "$b37" ]
			then
				$DIR/scripts/interactive_manh.py $curpeak_c.$curpeak_ps.peakdata.ld.annotated.assoc "$pvalcol" "$pscol" "$rscol" "$mafcol" "$chrcol"
			else
				$DIR/scripts/interactive_manh.b37.py $curpeak_c.$curpeak_ps.peakdata.ld.annotated.assoc "$pvalcol" "$pscol" "$rscol" "$mafcol" "$chrcol"
			fi
				
			assocfile=$mem

		else
			first=0
			echo "First peak at $curchr.$curps.$curs.$cura1.$cura2.$oldchr.$curpeak_p.$curpeak_rs.$curpeak_c.$curpeak_ps"
		fi
		oldpos=$curps
		oldchr=$curchr
		curpeak_p=$curpval
		curpeak_rs=$curs
		curpeak_c=$curchr
		curpeak_ps=$curps
		curpeak_a1=$cura1
		curpeak_a2=$cura2
		echo "Oldpos set to $curps $oldchr $curpeak_p $curpeak_rs $curpeak_c $curpeak_ps"
	else
		echo "Found another SNP, same peak"
		echo "Prev " $curpeak_rs $curpeak_ps $curpeak_c $curpeak_p
	# case where we have a high peak with multiple low pvals
	curpeak_rs=$(echo $line | awk -v rs=$rscoli -v curs=$curpeak_rs -v p=$pvalcoli -v cur=$curpeak_p '{if($p<cur){print $rs}else{print curs}}')
	curpeak_c=$(echo $line | awk -v p=$curpval -v cur=$curpeak_p -v chr=$chrcoli -v oldchr=$curpeak_c '{if(p<cur){print $chr}else{print oldchr}}')
	curpeak_ps=$(echo $line | awk -v p=$curpval -v cur=$curpeak_p -v ps=$pscoli -v oldps=$curpeak_ps '{if(p<cur){print $ps}else{print oldps}}')
	curpeak_a1=$(echo $line | awk -v p=$curpval -v cur=$curpeak_p -v ps=$a1coli -v oldps=$curpeak_a1 '{if(p<cur){print $ps}else{print oldps}}')
	curpeak_a2=$(echo $line | awk -v p=$curpval -v cur=$curpeak_p -v ps=$a2coli -v oldps=$curpeak_a2 '{if(p<cur){print $ps}else{print oldps}}')
	curpeak_p=$(echo $line | awk -v p=$pvalcoli -v cur=$curpeak_p '{if($p<cur){print $p}else{print cur}}')
		#echo "Post " $curpeak_rs $curpeak_ps $curpeak_c $curpeak_p

	#####
	## Here, we should also rename the file because as is, it remains the first signif seen.
	#####

	oldpos=$curpeak_ps
	oldchr=$curpeak_c
fi
echo $curpeak_rs > cprs
echo $curpeak_p > cpp
echo $curpeak_c > cpc
echo $curpeak_ps > cpps
echo $curpeak_a1 > cpa1
echo $curpeak_a2 > cpa2


done
curpeak_rs=$(cat cprs)
curpeak_p=$(cat cpp)
curpeak_c=$(cat cpc)
curpeak_ps=$(cat cpps)
curpeak_a1=$(cat cpa1)
curpeak_a2=$(cat cpa2)

if [ -z "$curpeak_rs" ]; then
	echo "No more peaks -DONE.["
	exit;
fi

echo -e "\n\n\nFINISHED READING FILE.\n\n\n"
echo curpeak_rs is \"$curpeak_rs\"

echo -e "Found peak at $(cat cprs) (p=$(cat cpp))"  
echo
echo -e "Fetching info from Ensembl..."
echo "chr"${curpeak_c}:${curpeak_ps}-${curpeak_ps}_${curpeak_a1}_${curpeak_a2}
			if [ -z "$b37" ]
			then
				$DIR/scripts/VarAnnot.b38.pl <(echo "chr"${curpeak_c}:${curpeak_ps}-${curpeak_ps}_${curpeak_a1}_${curpeak_a2}) > $curpeak_c.$curpeak_ps.line
			else
				$DIR/scripts/VarAnnot.pl <(echo "chr"${curpeak_c}:${curpeak_ps}-${curpeak_ps}_${curpeak_a1}_${curpeak_a2}) > $curpeak_c.$curpeak_ps.line
			fi
#$DIR/scripts/VarAnnot.b38.pl <(echo "chr"${curpeak_c}:${curpeak_ps}-${curpeak_ps}_${curpeak_a1}_${curpeak_a2}) > $curpeak_c.$curpeak_ps.line
#perl /nfs/team144/ds26/FunctionalAnnotation/20150505_new_tool/VarAnnot.pl <(echo "chr"${curpeak_c}:${curpeak_ps}-${curpeak_ps}_${curpeak_a1}_${curpeak_a2}) > $curpeak_c.$curpeak_ps.line
cp gwas_signals.tsv $curpeak_c.$curpeak_ps.signals
echo
echo -e "\t extracting $flank_bp" region...

cat <($cat $assocfile | head -n1 | tr '\t' ' ') <(awk -v chr=$chrcoli -v matchr=$curpeak_c -v ps=$pscoli -v matchps=$curpeak_ps -v flank=$flank_bp '$chr==matchr && $ps>(matchps - flank) && $ps<(matchps+flank)' <($cat $assocfile)) | sed 's/ /\t/g;s/\[b38\]//' > peakdata
awk -v cc=$chrcoli -v pc=$pscoli -v ic=$rscoli '{if(NR>1 && $ic!~/\:/ && $ic!~/rs/){$ic=$cc":"$pc}print}' peakdata > t
mv t peakdata






### UNCOMMENT THIS TO CLEANUP THE ASSOC FILE FOR LOW FREQ VAR HORIZ LINES



#/software/bin/Rscript --vanilla ~/association-scripts/cleanup.R
#			mv peakdata peakdata.unclean
#			mv peakdata.cleaned peakdata

cat  <(echo -e "snp\tchr\tpos") <(awk -v rs=$rscoli -v chr=$chrcoli -v ps=$pscoli 'BEGIN{OFS="\t"} NR>1 {print $rs, $chr, $ps}' peakdata) | tr ' ' '\t' > peakdata.chrpos
dbmeister.py --db $curpeak_c.$curpeak_ps.db --snp_pos peakdata.chrpos
dbmeister.py --db $curpeak_c.$curpeak_ps.db --refflat $REFFLAT
dbmeister.py --db $curpeak_c.$curpeak_ps.db --recomb_rate $RECOMB
for id in "${ids[@]}"
do 
			        sensible_start=$(($curpeak_ps-$flank_bp))
				if [ $sensible_start -le 0 ]
				then
				    sensible_start=1
				    echo -e "\n\n\nWARNING\t Negative start position changed to $sensible_start  : $curpeak_c $curpeak_ps (3)\n\n\n"
				fi

echo plink --bfile ${files[$id]} --chr $curpeak_c --from-bp $sensible_start --to-bp $(($curpeak_ps+$flank_bp)) --out $id --make-bed
plink --bfile ${files[$id]} --chr $curpeak_c --from-bp $sensible_start --to-bp $(($curpeak_ps+$flank_bp)) --out $id --make-bed
awk 'OFS="\t"{if($2!~/:/ && $2!~/rs/){$2="chr"$1":"$4}print}' $id.bim  | tr ';' '_'> t
mv t $id.bim

done
seq 0 $(($(echo $filelist | tr ',' '\n'| wc -l)-1)) > mergelist
plink --merge-list mergelist --make-bed --out merged --allow-no-sex
if [ -f merged-merge.missnp ]
	then
	grep -v -w -f merged-merge.missnp peakdata > t 
	mv t peakdata
	for id in "${ids[@]}"
	do 
	plink --bfile $id --exclude merged-merge.missnp --make-bed --out $id.tmp
	mv $id.tmp.bed $id.bed
	mv $id.tmp.bim $id.bim
	mv $id.tmp.fam $id.fam

done
plink --merge-list mergelist --make-bed --out merged --allow-no-sex
fi


awk '{if($1~/\:/ && !($1~/chr/)){$1="chr"$1}print}' peakdata.chrpos | tr '\t' ' '> t
mv t peakdata.chrpos
awk -v rspos=$rscoli '{if($rspos~/\:/ && !($rspos~/chr/)){$rspos="chr"$rspos}print}' peakdata | tr '\t' ' '> t
mv t peakdata
awk '{if($2~/\:/ && !($2~/chr/)){$2="chr"$2}print}' merged.bim |sed 's/\[b38\]//'> t
mv t merged.bim
refsnp=$(echo $curpeak_rs | awk '{if($0~/\:/ && !($0~/chr/)){$0="chr"$0}print}')
refsnp=$(echo $refsnp | awk -v c=$curpeak_c -v p=$curpeak_ps '{if($0!~/\:/ && $0!~/rs/){$0="chr"c":"p}print}')

echo 
echo
echo $curpeak_rs $refsnp
echo
echo
plink --bfile merged --r2 --ld-snp $refsnp  --ld-window-kb $ext_flank_kb $flank_kb_ext --ld-window 999999 --ld-window-r2 0 --out merged
cat <(echo -e "snp1\tsnp2\tdprime\trsquare") <(tail -n +2  merged.ld| awk 'OFS=" "{print $3,$6,$7,$7}') > t
mv t merged.ld
cp merged.ld $curpeak_c.$curpeak_ps.ld
echo
echo
echo $rscol
echo
echo
			        sensible_start=$(($curpeak_ps-$flank_bp))
				if [ $sensible_start -le 0 ]
				then
				    sensible_start=1
				    echo -e "\n\n\nWARNING\t Negative start position changed to $sensible_start  : $curpeak_c $curpeak_ps (4)\n\n\n"
				fi

echo locuszoom --metal peakdata --refsnp "$refsnp" --markercol "$rscol" --pvalcol "$pvalcol" --db $curpeak_c.$curpeak_ps.db --prefix $curpeak_c.$curpeak_ps.500kb --plotonly showAnnot=T showRefsnpAnnot=T annotPch="21,24,24,25,22,22,8,7" rfrows=20 geneFontSize=.4 --ld merged.ld --start=$sensible_start --end=$(($curpeak_ps+$flank_bp)) --chr=$curpeak_c showRecomb=T --delim ' '
locuszoom --metal peakdata --refsnp "$refsnp" --markercol "$rscol" --pvalcol "$pvalcol" --db $curpeak_c.$curpeak_ps.db --prefix $curpeak_c.$curpeak_ps.500kb --plotonly showAnnot=T showRefsnpAnnot=T annotPch="21,24,24,25,22,22,8,7" rfrows=20 geneFontSize=.4 --ld merged.ld --start=$sensible_start --end=$(($curpeak_ps+$flank_bp)) --chr=$curpeak_c showRecomb=T --delim ' '

#rm merged.* peakdata* 
#cat <(echo $(head -n1 peakdata) ld) <(join -1 $rscoli -2 1 <(sort -k$rscoli,$rscoli peakdata) <(cut -f2,3 -d' ' merged.ld |sort -k1,1) | perl -lane 'print join(@F[2..]), "\n";') > peakdata.ld
echo JOIN4
join --header -1 $rscoli -2 1 <(cat <(head -n1 peakdata) <(sort -k$rscoli,$rscoli peakdata)) <(cat <(echo $rscol ld) <(cut -f2,3 -d' ' merged.ld |sort -k1,1) ) > peakdata.ld
#join --header -j 1 <(cat <(head -n1 peakdata) <(sort -k$rscoli,$rscoli peakdata) | perl -lape 'print join($F[%ENV{"rscoli"}], @F), "\n";') <(cat <(echo $rscol ld) <(cut -f2,3 -d' ' merged.ld |sort -k1,1) ) > peakdata.ld

## Beware rsid is now number one and the columns are messed up.
cat=cat
assocfile=peakdata.ld
chrcoli=$(grep -w $chrcol <(paste <(seq 1 $($cat $assocfile | head -n1 | tr ' ' '\n'| wc -l)) <($cat $assocfile | head -n1 | tr ' ' '\n')) | cut -f1)
pscoli=$(grep -w $pscol <(paste <(seq 1 $($cat $assocfile | head -n1 | tr ' ' '\n'| wc -l)) <($cat $assocfile | head -n1 | tr ' ' '\n')) | cut -f1)
rscoli=$(grep -w $rscol <(paste <(seq 1 $($cat $assocfile | head -n1 | tr ' ' '\n'| wc -l)) <($cat $assocfile | head -n1 | tr ' ' '\n')) | cut -f1)
a1coli=$(grep -w $a1col <(paste <(seq 1 $($cat $assocfile | head -n1 | tr ' ' '\n'| wc -l)) <($cat $assocfile | head -n1 | tr ' ' '\n')) | cut -f1)
a2coli=$(grep -w $a2col <(paste <(seq 1 $($cat $assocfile | head -n1 | tr ' ' '\n'| wc -l)) <($cat $assocfile | head -n1 | tr ' ' '\n')) | cut -f1)
pvalcoli=$(grep -w $pvalcol <(paste <(seq 1 $($cat $assocfile | head -n1 | tr ' ' '\n'| wc -l)) <($cat $assocfile | head -n1 | tr ' ' '\n')) | cut -f1)
#cat <(echo $(head -n1 peakdata.ld) id gene consequence) <(join -1 3 -2 1 -a 1 <(sort -k3,3 peakdata.ld) <(~ag15/scripts/getrsid.pl <(awk '$NF>0.1' peakdata.ld | tail -n+2 | cut -f2-5 -d' ' )| sort -k1,1) | grep -v MARKER)| awk '{if(NF<22){$21="NA";$22="NA";$23="NA"}if(NR>1){un=$1;deux=$2;trois=$3;$1=deux;$2=trois;$3=un;}print}' > peakdata.ld.annotated
expnumcol=$(($(head -n1 peakdata.ld | tr ' ' '\n' | wc -l)+3))


if [ -z "$b37" ]
then
	join -a 1 -1 $pscoli -2 1 --header <(cat <(head -n1 peakdata.ld) <(tail -n+2 peakdata.ld | sort -k$pscoli,$pscoli)) <( cat <(echo $pscol id gene consequence) <($DIR/scripts/getrsid.b38.pl <(awk '$NF>0.1' peakdata.ld | tail -n+2 | awk -v chr=$chrcoli -v ps=$pscoli -v a1=$a1coli -v a2=$a2coli '{print $chr,$ps,$a1,$a2}') | sort -k1,1))  |awk -v en=$expnumcol '{if(NF<en){$(en-2)="NA";$(en-1)="NA";$en="NA"}print}'> peakdata.ld.annotated
	$DIR/scripts/get_phenotype.b38.pl <(cut -f$(($expnumcol-2)) -d' ' peakdata.ld.annotated | grep -v -e NA -e novel -e id) | sed 's/ /&/g'| sed 's/\t/\t"/;s/$/"/' > phenotypes
else
	join -a 1 -1 $pscoli -2 1 --header <(cat <(head -n1 peakdata.ld) <(tail -n+2 peakdata.ld | sort -k$pscoli,$pscoli)) <( cat <(echo $pscol id gene consequence) <($DIR/scripts/getrsid.b38.pl b37 <(awk '$NF>0.1' peakdata.ld | tail -n+2 | awk -v chr=$chrcoli -v ps=$pscoli -v a1=$a1coli -v a2=$a2coli '{print $chr,$ps,$a1,$a2}') | sort -k1,1))  |awk -v en=$expnumcol '{if(NF<en){$(en-2)="NA";$(en-1)="NA";$en="NA"}print}'> peakdata.ld.annotated
	$DIR/scripts/get_phenotype.b38.pl b37 <(cut -f$(($expnumcol-2)) -d' ' peakdata.ld.annotated | grep -v -e NA -e novel -e id) | sed 's/ /&/g'| sed 's/\t/\t"/;s/$/"/' > phenotypes
fi

#cat <(echo $(head -n1 peakdata.ld.annotated) assoc) <(join -1 21 -2 1 -a 1 <(tail -n+2 peakdata.ld.annotated | sort -k21,21) <(sort -k1,1 phenotypes) |awk '{if(NF!=24){$24="NA"}print}'| perl -lane 'print join(" ", @F[1..20], $F[0],@F[21..24]);') |sed 's/ /,/g;s/_/ /g'> $curpeak_c.$curpeak_ps.peakdata.ld.annotated.assoc
echo JOIN6
join --header -a 1 -1 $(($expnumcol-2)) -2 1 <(cat <(head -n1 peakdata.ld.annotated) <(tail -n+2 peakdata.ld.annotated|  sort -k$(($expnumcol-2)),$(($expnumcol-2)))) <( cat <(echo id assoc) <(sort -k1,1 phenotypes)) |awk -v field=$((expnumcol+2)) '{if(NF!=field){$field="NA"}print}' |sed 's/ /,/g;s/&/ /g'> $curpeak_c.$curpeak_ps.peakdata.ld.annotated.assoc

if [ -z "$b37" ]
then
	$DIR/scripts/interactive_manh.py $curpeak_c.$curpeak_ps.peakdata.ld.annotated.assoc "$pvalcol" "$pscol" "$rscol" "$mafcol" "$chrcol"
else
	$DIR/scripts/interactive_manh.b37.py $curpeak_c.$curpeak_ps.peakdata.ld.annotated.assoc "$pvalcol" "$pscol" "$rscol" "$mafcol" "$chrcol"
fi

for i in `ls *.bak`
do
	echo "Back-plotting file $i."
	snp=$(echo $i | sed 's/peakdata.//;s/.bak//')
	#chr=$(echo $snp | sed 's/chr//;s/\:.*//')
	#chr=$(cut -f$chrcoli -d ' ' $i | sort -u | grep -vw 'chr' | head -1 | sed 's/chr//')
	mem=$assocfile
	assocfile=$i
			chrcoli=$(grep -w $chrcol <(paste <(seq 1 $(cat $assocfile | head -n1 | tr ' ' '\n'| wc -l)) <(cat $assocfile | head -n1 | tr ' ' '\n')) | cut -f1)
			pscoli=$(grep -w $pscol <(paste <(seq 1 $(cat $assocfile | head -n1 | tr ' ' '\n'| wc -l)) <(cat $assocfile | head -n1 | tr ' ' '\n')) | cut -f1)
	assocfile=$mem
	pos=$(echo $snp | sed 's/.*\://')
	posa=$(grep $pos $i | cut -f$pscoli -d' ')
	chr=$(grep $pos $i | cut -f$chrcoli -d' ')
	echo "snp $snp chr $chr pos $posa"
	pos=$posa
	sensible_start=$(($pos-$flank_bp))
	if [ $sensible_start -le 0 ]
	then
	    sensible_start=1
	    echo -e "\n\n\nWARNING\t Negative start position changed to $sensible_start  : $curpeak_c $curpeak_ps $pos $snp $chr (5)\n\n\n"
	fi
	if [ -z "$b37" ]
	then
		$DIR/scripts/interactive_manh.py $chr.$pos.peakdata.ld.annotated.assoc "$pvalcol" "$pscol" "$rscol" "$mafcol" "$chrcol"
	else
		$DIR/scripts/interactive_manh.b37.py $chr.$pos.peakdata.ld.annotated.assoc "$pvalcol" "$pscol" "$rscol" "$mafcol" "$chrcol"
	fi
	/nfs/team144/software/locuszoom-1.3/locuszoom/bin/locuszoom --metal $i --refsnp "$snp" --markercol "$rscol" --pvalcol "$pvalcol" --db $chr.$pos.db --prefix $chr.$pos.500kb --plotonly showAnnot=T showRefsnpAnnot=T annotPch="21,24,24,25,22,22,8,7" rfrows=20 geneFontSize=.4 --ld $chr.$pos.ld --start=$sensible_start --end=$(($pos+$flank_bp)) --chr=$chr showRecomb=T --delim ' '

done
rm cp* merge* peak* 0.* *.db *signal* *.line *.ld 

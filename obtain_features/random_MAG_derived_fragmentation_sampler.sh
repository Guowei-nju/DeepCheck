# This is used on a genome where random fragmentation has been applied to generate varying levels of completeness/contamination
# first argument is genome file
# get base pair number
round() {
  printf "%.${2}f" "${1}"
}
charcount=`grep -v ">" $1 | wc | awk '{print $3-$1}'`
echo "Base pair count is: " $charcount "(" $((charcount/1000)) " Mbp )"
totalNumberOfContigs=`grep ">" $1 | wc -l`
onePercentOfContigs=`bc -l <<< "scale=4; ($totalNumberOfContigs/100)"` 
onePercentOfContigsINT=$(round ${onePercentOfContigs} 0)

echo "OnePC is: " $onePercentOfContigsINT
tenPercentOfContigs=`bc -l <<< "scale=4; ($totalNumberOfContigs/10)"` 
tenPercentOfContigsINT=${tenPercentOfContigs%%.*}
genomeName=$1
#remove everything up to and including  last / (ie take out path)
genomeName=${genomeName##*/}
# remove .fna extention
genomeName=${genomeName::-4}
echo "Calling genome "$genomeName
#genome_chunk=$((charcount/chunkcount))
# echo "Splitting genome files into chunks of" $genome_chunk "base-pairs..."
echo "Randomly sampling genome at 5% completeness intervals pieces"
#module load bbmap
# `reformat.sh in=$1 out=./fragged_genome.fna breaklength=$genome_chunk` >> log.txt xxxxrxe
rm -rf "fragged_"$genomeName
mkdir "fragged_"$genomeName
cd "fragged_"$genomeName
cat $1 > $genomeName".fna"
index=0
samplerate=100
multiplier=10
while [ $index -le 15 ] # only do down to 15 % completeness in 5 % increments
do
  echo "Processing genome #" $index 
  #make new genomes from random sampling
  #samplerate=`bc -l <<< "$samplerate / 100"` sampleseed=1 ignorejunk=True    $(round ${$((onePercentOfContigsINT*samplerate))} 0)
  target=`bc -l <<< "scale=4; ($onePercentOfContigs*$samplerate)"` 
  echo 'TAGET IS: ' $target
  reformat.sh in=$genomeName".fna" out=$genomeName"_"$samplerate"_pc_complete.fna" samplereadstarget=$(round ${target} 0)
  # now check the actual completeness percentage
  charcount_fragged=`grep -v ">" $genomeName"_"$samplerate"_pc_complete.fna" | wc | awk '{print $3-$1}'`
  echo "Charcount_fragged: " $charcount_fragged "charcount: " $charcount 
  actualCompPercentage=`bc -l <<< "scale=6; (($charcount_fragged / $charcount) * 100)"` # we use this to label
  cat $genomeName"_"$samplerate"_pc_complete.fna" > TMP_.fna
  echo "Actual Comp is" $actualCompPercentage
  rm $genomeName"_"$samplerate"_pc_complete.fna"
  cat TMP_.fna > $genomeName"_"$actualCompPercentage"_pc_complete.fna"
  rm TMP_.fna
  

  #generate contamination
  if (($samplerate >  30)); then
    #make contamination source from incomplete genome:
    cat $genomeName"_"$actualCompPercentage"_pc_complete.fna" > cont_source.fna
    f2pc=`bc -l <<< "scale=4; ($onePercentOfContigs * 2)"`
    f5pc=`bc -l <<< "scale=4; ($onePercentOfContigs * 5)"`
    f8pc=`bc -l <<< "scale=4; ($onePercentOfContigs * 8)"`
    f10pc=`bc -l <<< "scale=4; ($onePercentOfContigs * 10)"`
    f15pc=`bc -l <<< "scale=4; ($onePercentOfContigs * 15)"`
    f20pc=`bc -l <<< "scale=4; ($onePercentOfContigs * 20)"`
    f25pc=`bc -l <<< "scale=4; ($onePercentOfContigs * 25)"`
    f30pc=`bc -l <<< "scale=4; ($onePercentOfContigs * 30)"`
    int2pc=${f2pc%%.*}
    int5pc=${f5pc%%.*}
    int8pc=${f8pc%%.*}
    int10pc=${f10pc%%.*}
    int15pc=${f15pc%%.*}
    int20pc=${f20pc%%.*}
    int25pc=${f25pc%%.*}
    int30pc=${f30pc%%.*}
    reformat.sh in=cont_source.fna out=2pc.fna samplereadstarget=$int2pc ignorejunk=True
    reformat.sh in=cont_source.fna out=5pc.fna samplereadstarget=$int5pc ignorejunk=True
    reformat.sh in=cont_source.fna out=8pc.fna samplereadstarget=$int8pc ignorejunk=True
    reformat.sh in=cont_source.fna out=10pc.fna samplereadstarget=$int10pc ignorejunk=True
    reformat.sh in=cont_source.fna out=15pc.fna samplereadstarget=$int15pc ignorejunk=True
    reformat.sh in=cont_source.fna out=20pc.fna samplereadstarget=$int20pc ignorejunk=True
    reformat.sh in=cont_source.fna out=25pc.fna samplereadstarget=$int25pc ignorejunk=True
    reformat.sh in=cont_source.fna out=30pc.fna samplereadstarget=$int30pc ignorejunk=True
    rm cont_source.fna

    # now calculate the size of contamination chunks
    f2pcCont=`grep -v ">" 2pc.fna | wc | awk '{print $3-$1}'`
    f2pcContActual=`bc -l <<< "scale=6; (($f2pcCont / $charcount) * 100)"` # actual basepairs of '2 percent contamination'
    
    f5pcCont=`grep -v ">" 5pc.fna | wc | awk '{print $3-$1}'`
    f5pcContActual=`bc -l <<< "scale=6; (($f5pcCont / $charcount) * 100)"` # actual basepairs of '5 percent contamination'
    
    f8pcCont=`grep -v ">" 8pc.fna | wc | awk '{print $3-$1}'`
    f8pcContActual=`bc -l <<< "scale=6; (($f8pcCont / $charcount) * 100)"` # actual basepairs of '8 percent contamination'
    
    f10pcCont=`grep -v ">" 10pc.fna | wc | awk '{print $3-$1}'`
    f10pcContActual=`bc -l <<< "scale=6; (($f10pcCont / $charcount) * 100)"` # actual basepairs of '10 percent contamination'
    
    f15pcCont=`grep -v ">" 15pc.fna | wc | awk '{print $3-$1}'`
    f15pcContActual=`bc -l <<< "scale=6; (($f15pcCont / $charcount) * 100)"` # actual basepairs of '15 percent contamination'
    
    f20pcCont=`grep -v ">" 20pc.fna | wc | awk '{print $3-$1}'`
    f20pcContActual=`bc -l <<< "scale=6; (($f20pcCont / $charcount) * 100)"` # actual basepairs of '20 percent contamination'
    
    f25pcCont=`grep -v ">" 25pc.fna | wc | awk '{print $3-$1}'`
    f25pcContActual=`bc -l <<< "scale=6; (($f25pcCont / $charcount) * 100)"` # actual basepairs of '20 percent contamination'
    
    f30pcCont=`grep -v ">" 30pc.fna | wc | awk '{print $3-$1}'`
    f30pcContActual=`bc -l <<< "scale=6; (($f30pcCont / $charcount) * 100)"` # actual basepairs of '20 percent contamination'

    
    cat $genomeName"_"$actualCompPercentage"_pc_complete.fna" 2pc.fna  > $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f2pcContActual"pc_contaminated.fna"
    cat $genomeName"_"$actualCompPercentage"_pc_complete.fna" 5pc.fna  > $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f5pcContActual"pc_contaminated.fna"
    cat $genomeName"_"$actualCompPercentage"_pc_complete.fna" 8pc.fna  > $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f8pcContActual"pc_contaminated.fna"
    cat $genomeName"_"$actualCompPercentage"_pc_complete.fna" 10pc.fna  > $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f10pcContActual"pc_contaminated.fna"
    cat $genomeName"_"$actualCompPercentage"_pc_complete.fna" 15pc.fna  > $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f15pcContActual"pc_contaminated.fna"
    cat $genomeName"_"$actualCompPercentage"_pc_complete.fna" 20pc.fna  > $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f20pcContActual"pc_contaminated.fna"
    cat $genomeName"_"$actualCompPercentage"_pc_complete.fna" 25pc.fna  > $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f25pcContActual"pc_contaminated.fna"
    cat $genomeName"_"$actualCompPercentage"_pc_complete.fna" 30pc.fna  > $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f30pcContActual"pc_contaminated.fna"
    rm 2pc.fna
    rm 5pc.fna
    rm 8pc.fna
    rm 10pc.fna
    rm 15pc.fna
    rm 20pc.fna
    rm 25pc.fna
    rm 30pc.fna
    rm cont_source.fna
  fi
  cat $genomeName"_"$actualCompPercentage"_pc_complete.fna" > $genomeName"_"$actualCompPercentage"_pc_completeXXX0pc_contaminated.fna"
  rm $genomeName"_"$actualCompPercentage"_pc_complete.fna"
  rm shredded_whole.fna
  if (($samplerate <  50)); then
    rm $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f30pcContActual"pc_contaminated.fna"
    rm $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f25pcContActual"pc_contaminated.fna"
  fi
  if (($samplerate <  40)); then
    rm $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f10pcContActual"pc_contaminated.faa"
    rm $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f8pcContActual"pc_contaminated.faa"
  fi
  index=$((index+1))
  multiplier=$((multiplier-1))
  if (($samplerate > 80)); then
    samplerate=$((samplerate-2))
  else
    samplerate=$((samplerate-5))
  fi
done
rm $genomeName".fna"
#for i in *.fna; do sed -e '/>/s/^/@/' -e '/>/s/$/#/' $i | tr -d "\n" | tr "@" "\n" | sort -t "_" -k2n | tr "#" "\n" | sed -e '/^$/d'^C > tmp && mv tmp $i; done
#ls *.fna > fnalist.txt
#for i in *.fna; do sed -e "s/>/>$i/g" -i $i; awk '/>/{print $0(++i)}!/>/' $i > tmp && mv tmp $i; done
#module load fastani
#echo "Calling FastANI on protein level."
#`fastANI --ql fnalist.txt --rl fnalist.txt --fragLen 200 -o fastANI_fna_output.tsv &>> log.txt`
#echo "Calling FastANI on nucleotide level."
#`fastANI --ql fnalist.txt --rl fnalist.txt --fragLen 200 -o fastANI_FNA_output.tsv &>> log.txt`
#module unload fastANI
for i in *.fna; do sed '/^>/ s/ .*//' -i $i; done
for i in *.fna; do awk '/^>/{$0=$0"_"(++i)}1' $i > tmp"$i" && mv tmp"$i" $i; done
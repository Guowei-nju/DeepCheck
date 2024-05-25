# first argument is genome file
# get base pair number
charcount=`grep -v ">" $1 | wc | awk '{print $3-$1}'`
numer=$((charcount+10000))
totalContigsWhen20kbSheared=$((numer/19400))
#this is for contamination: 
chunkcount=100
genome_chunk=$((charcount/chunkcount))
echo "Base pair count is: " $charcount "(" $((charcount/1000)) " Mbp )"
# get name
genomeName=$1
#remove everything up to and including  last /
genomeName=${genomeName##*/}
# remove .fna extention
genomeName=${genomeName::-4}
echo "Calling genome "$genomeName
onePercentNumer=$((charcount+50))
onePercentOfContigs=$((totalContigsWhen20kbSheared+$onePercentNumer/100))
#this is for cont
chunkcount=100
#genome_chunk=$((charcount/chunkcount))
# echo "Splitting genome files into chunks of" $genome_chunk "base-pairs..."
echo "Randomly sampling genome at 5% completeness intervals using 20kb pieces"
#module load bbmap
# `reformat.sh in=$1 out=./fragged_genome.fna breaklength=$genome_chunk` >> log.txt
rm -rf "fragged_"$genomeName
mkdir "fragged_"$genomeName
cd "fragged_"$genomeName
cat $1 > ./$genomeName".fna"
# let's deal with any tiny snippets of DNA by incorporating (mostly plasmid) dna into core genome. Not the best way perhaps but at least the most consistent. 
#grep -v ">" $genomeName".fna" > $genomeName".fna_temp"
#sed  -i "1i >$i" $genomeName".fna_temp"
#cat $genomeName".fna_temp" > $genomeName".fna"
#rm $genomeName".fna_temp"
index=0
samplerate=100
while [ $index -le 10 ] # go down in 5 % increments
do
  echo "Processing genome #" $index
  shred.sh in=$genomeName".fna" out=shredded_whole.fna median=20000 minlength=2500 variance=2000
  echo "nihao #" $index
  # this is an exact number of contigs needed to get to the samplerate of interest
  currentchunk=`bc -l <<< "scale=4; (($totalContigsWhen20kbSheared/100) * $samplerate)"`
  # this is how many contigs we need for a given samplerate, rounded to int
  chunkINT=${currentchunk%%.*} # we use this to sample
  # now let's determine what pc comp we've generated using the closest whole contig number: 
  #actualCompPercentage0=`bc -l <<< "scale=6; ((($chunkINT * 20000)/ $charcount) * 100)"` # we use this to label
  #shred.sh in=$genomeName".fna" out=shredded_whole.fna length=20000 minlength=2500
  shuffle.sh in=shredded_whole.fna out=shredded_whole2.fna
  # it is now shuffled for this round of sampling 
  #reformat.sh in=$genomeName".fna" out=$genomeName"_"$samplerate"_pc_complete.fna" forcetrimleft=$((fivepercent_chunk*index)) ?// this is for trimming
  reformat.sh in=shredded_whole2.fna out=$genomeName"_"$samplerate"_pc_complete.fna" samplereadstarget=$chunkINT ignorejunk=t
  charcount_fragged=`grep -v ">" $genomeName"_"$samplerate"_pc_complete.fna" | wc | awk '{print $3-$1}'`
  echo "Charcount_fragged: " $charcount_fragged "charcount: " $charcount
  actualCompPercentage=`bc -l <<< "scale=6; (($charcount_fragged / $charcount) * 100)"` # we use this to label
  cat $genomeName"_"$samplerate"_pc_complete.fna" > TMP_.fna
  echo "Actual Comp is" $actualCompPercentage
  rm $genomeName"_"$samplerate"_pc_complete.fna"
  cat TMP_.fna > $genomeName"_"$actualCompPercentage"_pc_complete.fna"
  rm TMP_.fna
  rm shredded_whole2.fna
  
  # at this point it's worth checking if we've generated a 100% complete (fragged) genome (or at least 99%), and if we haven't, make one. 
  # together with some contamination! 
  # also we will import the unfragged genomes from the trimmed cycle, but that's a somewhat different story 
  if (($samplerate == 100)); then 
      cat $genomeName".fna" > cont_source.fna
      shred.sh in=cont_source.fna out=shredded.fna median=20000 minlength=2500 variance=5 #length=$((genome_chunk+1))
      f2pc=`bc -l <<< "scale=4; (($totalContigsWhen20kbSheared/100) * 2)"`
      f5pc=`bc -l <<< "scale=4; (($totalContigsWhen20kbSheared/100) * 5)"`
      f8pc=`bc -l <<< "scale=4; (($totalContigsWhen20kbSheared/100) * 8)"`
      f10pc=`bc -l <<< "scale=4; (($totalContigsWhen20kbSheared/100) * 10)"`
      f15pc=`bc -l <<< "scale=4; (($totalContigsWhen20kbSheared/100) * 15)"`
      f20pc=`bc -l <<< "scale=4; (($totalContigsWhen20kbSheared/100) * 20)"`
      f25pc=`bc -l <<< "scale=4; (($totalContigsWhen20kbSheared/100) * 25)"`
      f30pc=`bc -l <<< "scale=4; (($totalContigsWhen20kbSheared/100) * 30)"`
      int2pc=${f2pc%%.*}
      int5pc=${f5pc%%.*}
      int8pc=${f8pc%%.*}
      int10pc=${f10pc%%.*}
      int15pc=${f15pc%%.*}
      int20pc=${f20pc%%.*}
      int25pc=${f25pc%%.*}
      int30pc=${f30pc%%.*}
      reformat.sh in=shredded.fna out=2pc.fna samplereadstarget=$int2pc ignorejunk=t
      shuffle.sh in=shredded.fna out=shredded2.fna 
      reformat.sh in=shredded2.fna out=5pc.fna samplereadstarget=$int5pc ignorejunk=t
      shuffle.sh in=shredded2.fna out=shredded3.fna 
      reformat.sh in=shredded3.fna out=8pc.fna samplereadstarget=$int8pc ignorejunk=t
      shuffle.sh in=shredded3.fna out=shredded4.fna 
      reformat.sh in=shredded4.fna out=10pc.fna samplereadstarget=$int10pc ignorejunk=t
      shuffle.sh in=shredded4.fna out=shredded5.fna 
      reformat.sh in=shredded5.fna out=15pc.fna samplereadstarget=$int15pc ignorejunk=t
    
      shuffle.sh in=shredded5.fna out=shredded6.fna 
      reformat.sh in=shredded6.fna out=20pc.fna samplereadstarget=$int20pc ignorejunk=t
    
      shuffle.sh in=shredded6.fna out=shredded7.fna 
      reformat.sh in=shredded7.fna out=25pc.fna samplereadstarget=$int25pc ignorejunk=t
    
      shuffle.sh in=shredded7.fna out=shredded8.fna 
      reformat.sh in=shredded8.fna out=30pc.fna samplereadstarget=$int30pc ignorejunk=t
      rm shredded.fna
      rm shredded2.fna
      rm shredded3.fna
      rm shredded4.fna
      rm shredded5.fna
      rm shredded6.fna
      rm shredded7.fna
      rm shredded8.fna
      
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
      f25pcContActual=`bc -l <<< "scale=6; (($f25pcCont / $charcount) * 100)"` # actual basepairs of '25 percent contamination'
    
      f30pcCont=`grep -v ">" 30pc.fna | wc | awk '{print $3-$1}'`
      f30pcContActual=`bc -l <<< "scale=6; (($f30pcCont / $charcount) * 100)"` # actual basepairs of '30 percent contamination'
      
      cat cont_source.fna 2pc.fna  > $genomeName"_100_pc_completeXXX"$f2pcContActual"pc_contaminated.fna"
      cat cont_source.fna 5pc.fna  > $genomeName"_100_pc_completeXXX"$f5pcContActual"pc_contaminated.fna"
      cat cont_source.fna 8pc.fna  > $genomeName"_100_pc_completeXXX"$f8pcContActual"pc_contaminated.fna"
      cat cont_source.fna 10pc.fna  > $genomeName"_100_pc_completeXXX"$f10pcContActual"pc_contaminated.fna"
      cat cont_source.fna 15pc.fna  > $genomeName"_100_pc_completeXXX"$f15pcContActual"pc_contaminated.fna"
      cat cont_source.fna 20pc.fna  > $genomeName"_100_pc_completeXXX"$f20pcContActual"pc_contaminated.fna"
      cat cont_source.fna 25pc.fna  > $genomeName"_100_pc_completeXXX"$f25pcContActual"pc_contaminated.fna"
      cat cont_source.fna 30pc.fna  > $genomeName"_100_pc_completeXXX"$f30pcContActual"pc_contaminated.fna"
      cat cont_source.fna > $genomeName"_100_pc_completeXXX0pc_contaminated.fna"
      rm 2pc.fna
      rm 5pc.fna
      rm 8pc.fna
      rm 10pc.fna
      rm 15pc.fna
      rm 20pc.fna
      rm 30pc.fna
      rm 25pc.fna
      rm cont_source.fna
      
    #cat $genomeName"_"$actualCompPercentage"_pc_complete.fna" > $genomeName"_"$actualCompPercentage"_pc_completeXXX0pc_contaminated.fna"
    rm $genomeName"_"$actualCompPercentage"_pc_complete.fna"
    rm shredded_whole.fna
  fi
  
  if (($samplerate <  100)); then
    if (($samplerate >  30)); then
      #make contamination source
      grep -v ">" $genomeName"_"$actualCompPercentage"_pc_complete.fna" > cont_source.fna
      sed  -i "1i >$i" cont_source.fna
      shred.sh in=cont_source.fna out=shredded.fna median=20000 minlength=2500 variance=2000 #length=$((genome_chunk+1))
      f2pc=`bc -l <<< "scale=4; (($totalContigsWhen20kbSheared/100) * 2)"`
      f5pc=`bc -l <<< "scale=4; (($totalContigsWhen20kbSheared/100) * 5)"`
      f8pc=`bc -l <<< "scale=4; (($totalContigsWhen20kbSheared/100) * 8)"`
      f10pc=`bc -l <<< "scale=4; (($totalContigsWhen20kbSheared/100) * 10)"`
      f15pc=`bc -l <<< "scale=4; (($totalContigsWhen20kbSheared/100) * 15)"`
      f20pc=`bc -l <<< "scale=4; (($totalContigsWhen20kbSheared/100) * 20)"`
      f25pc=`bc -l <<< "scale=4; (($totalContigsWhen20kbSheared/100) * 25)"`
      f30pc=`bc -l <<< "scale=4; (($totalContigsWhen20kbSheared/100) * 30)"`
      int2pc=${f2pc%%.*}
      int5pc=${f5pc%%.*}
      int8pc=${f8pc%%.*}
      int10pc=${f10pc%%.*}
      int15pc=${f15pc%%.*}
      int20pc=${f20pc%%.*}
      int25pc=${f25pc%%.*}
      int30pc=${f30pc%%.*}
    
      reformat.sh in=shredded.fna out=2pc.fna samplereadstarget=$int2pc ignorejunk=t
      shuffle.sh in=shredded.fna out=shredded2.fna 
      reformat.sh in=shredded2.fna out=5pc.fna samplereadstarget=$int5pc ignorejunk=t
      shuffle.sh in=shredded2.fna out=shredded3.fna 
      reformat.sh in=shredded3.fna out=8pc.fna samplereadstarget=$int8pc ignorejunk=t
      shuffle.sh in=shredded3.fna out=shredded4.fna 
      reformat.sh in=shredded4.fna out=10pc.fna samplereadstarget=$int10pc ignorejunk=t
      shuffle.sh in=shredded4.fna out=shredded5.fna 
      reformat.sh in=shredded5.fna out=15pc.fna samplereadstarget=$int15pc ignorejunk=t
    
      shuffle.sh in=shredded5.fna out=shredded6.fna 
      reformat.sh in=shredded6.fna out=20pc.fna samplereadstarget=$int20pc ignorejunk=t
    
      shuffle.sh in=shredded6.fna out=shredded7.fna 
      reformat.sh in=shredded7.fna out=25pc.fna samplereadstarget=$int25pc ignorejunk=t
    
      shuffle.sh in=shredded7.fna out=shredded8.fna 
      reformat.sh in=shredded8.fna out=30pc.fna samplereadstarget=$int30pc ignorejunk=t
      rm shredded.fna
      rm shredded2.fna
      rm shredded3.fna
      rm shredded4.fna
      rm shredded5.fna
      rm shredded6.fna
      rm shredded7.fna
      rm shredded8.fna
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
      f25pcContActual=`bc -l <<< "scale=6; (($f25pcCont / $charcount) * 100)"` # actual basepairs of '25 percent contamination'
    
      f30pcCont=`grep -v ">" 30pc.fna | wc | awk '{print $3-$1}'`
      f30pcContActual=`bc -l <<< "scale=6; (($f30pcCont / $charcount) * 100)"` # actual basepairs of '30 percent contamination'
    
    
    
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
    
      if (($samplerate <  40)); then
        rm $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f30pcContActual"pc_contaminated.fna"
        rm $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f25pcContActual"pc_contaminated.fna"
      fi
    fi
  
    cat $genomeName"_"$actualCompPercentage"_pc_complete.fna" > $genomeName"_"$actualCompPercentage"_pc_completeXXX0pc_contaminated.fna"
    rm $genomeName"_"$actualCompPercentage"_pc_complete.fna"
    rm shredded_whole.fna
  fi
  index=$((index+1))
  samplerate=$((samplerate-5))
done
rm $genomeName".fna"
for i in *.fna; do sed '/^>/ s/ .*//' -i $i & done
for i in *.fna; do awk '/^>/{$0=$0"_"(++i)}1' $i > tmp"$i" && mv tmp"$i" $i; done
mkdir ../all_genomes
cp *.fna ../all_genomes

#!/bin/bash

# this script will take a reference (or assume human) and map the graph as well as final assembly to it
# assumes HPC version of reference is named same but .hpc.fasta, if it doesn't exist it will be created
# the output is a csv file suitable for bandage to new chromosome assignments
# it will also output assignments of the final sequences to each chromosome 

if [ x$SLURM_CPUS_PER_TASK = x -o x$SLURM_CPUS_PER_TASK = xundefined -o x$SLURM_CPUS_PER_TASK = x0 ]; then
  cores=`grep -c ^processor /proc/cpuinfo`
else
  cores=$SLURM_CPUS_PER_TASK
fi

p=`pwd`

mashmap=$(which mashmap 2>/dev/null)
if [ "x$mashmap" == "x" ]; then
   mashmap="/data/korens/devel/MashMap/build/bin/mashmap"
fi
if [ "x$mashmap" == "x" ]; then
   echo "Error: mashmap not found"
   exit 1
fi 

seqtk=$(which seqtk 2>/dev/null)
if [ "x$seqtk" == "x" ]; then
   module load seqtk
   seqtk=$(which seqtk 2>/dev/null)
fi
if [ "x$seqtk" == "x" ]; then
   echo "Error: seqtk not found"
   exit
fi

neigh=$(which neighborhood 2>/dev/null)
if [ "x$neigh" == "x" ]; then
   neigh="/gfacpp/build/neighborhood"
fi
if [ "x$neigh" == "x" ]; then
   echo "Error: gfacpp not found"
   exit
fi

ref="/data/Phillippy/t2t-share/assemblies/release/v2.0/chm13v2.0.fasta"
if [[ "$#" -ge 1 ]]; then
   echo "Using custom reference sequence $1"
   ref=`realpath $1`
   hpcRef=`echo "$ref"|sed s/.fasta/.hpc.fasta/g`
   if [ ! -e $hpcRef ]; then
      $seqtk hpc $ref > $hpcRef
   fi
fi
hpcRef=`echo "$ref" |sed s/.fasta/.hpc.fasta/g`

if [ ! -e assembly.mashmap.out ]; then
   $mashmap -r $ref -q assembly.fasta --pi 95 -s 10000 -t 8 -o assembly.mashmap.out
fi

if [ ! -e compressed.mashmap.out ]; then
   isRUK=`ls 6-rukki/unitig*.fasta 2>/dev/null |wc -l |awk '{print $1}'`
   isHIC=`ls 8-hicPipeline/unitigs.hpc.fasta 2>/dev/null |wc -l |awk '{print $1}'`

   if [ $isRUK -gt 0 ]; then
      cd 6-rukki
      $mashmap -r $hpcRef -q unitig*.fasta --pi 95 -s 10000 -t 8
      cd ..
      ln -s 6-rukki/mashmap.out compressed.mashmap.out
   elif [ $isHIC -gt 0 ]; then
      cd 8-hicPipeline
      $mashmap -r $hpcRef -q unitigs.hpc.fasta --pi 95 -s 10000 -t 8
      cd ..
      ln -s 8-hicPipeline/mashmap.out compressed.mashmap.out
   else
      cd 5-untip
	  if [ ! -e unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.fasta ]; then 
	  	  cat unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.gfa |awk '{if (match($1, "^S")) { print ">"$2; print $3}}'|fold -c  > unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.fasta
	  fi
	  $mashmap -r $hpcRef -q unitig*.fasta --pi 95 -s 10000 -t 8
	  cd ..
	  ln -s 5-untip/mashmap.out compressed.mashmap.out
   fi
fi

hasNC=`cat assembly.mashmap.out |grep -c NC`
g="."
if [ $hasNC -gt 0 ]; then
   g="NC"
fi
hasCNC=`cat compressed.mashmap.out |grep -c NC`
cg="."
if [ $hasCNC -gt 0 ]; then
   cg="NC"
fi
label1="mat-"
label2="pat-"

isMAT=`grep -c "mat-" assembly.fasta`
if [ $isMAT -eq 0 ]; then
   isHF=`grep -c "h1tg" assembly.fasta`
   if [ $isHF -eq 0 ]; then
      isUnpaired=`grep -c "contig-" assembly.fasta`
      if [ $isUnpaired -eq 0 ]; then
         label1="haplotype1-"
         label2="haplotype2-"
      else
         label1="contig-"
         label2="none-ignore"
      fi
   else
      label1="h1tg"
      label2="h2tg"
   fi
fi
echo "$isMAT $label1 $label2 compNC: $cg regNC: $g"

minLen=5000000
echo -e "node\tchr" > assembly.homopolymer-compressed.chr.csv
for i in `cat compressed.mashmap.out |awk '{if ($NF > 99) print $6}'|sort |uniq`; do
   echo "Chr $i"
   cat compressed.mashmap.out |awk -v M=$minLen '{if ($NF > 99 && $4-$3 > M) print $1"\t"$6"\t"$2}'|grep -w $i |sort -srnk3,3 |awk '{print $1"\t"$2}' |sort |uniq | grep "$cg" >> assembly.homopolymer-compressed.chr.csv
done

cat assembly.homopolymer-compressed.chr.csv |awk '{print $1}' |sort |uniq > tmp4
$neigh assembly.homopolymer-compressed.noseq.gfa tmp.gfa -n tmp4 --drop-sequence -r 1000
cat tmp.gfa |grep "^S" |awk '{print $2}' > tmp4

#second pass to get missing chr w/shorter matches
minLen=`echo $minLen |awk '{print $1/10}'`
cat compressed.mashmap.out |grep -w -v -f tmp4 |awk -v M=$minLen '{if ($NF > 99 && $4-$3 > M) print $1"\t"$6}' |sort |uniq | grep "$cg" >> assembly.homopolymer-compressed.chr.csv
rm tmp.gfa tmp4

minLen=`echo $minLen |awk '{print $1*10/5}'`
cat assembly.mashmap.out |grep "$label1" |grep $g |awk -v M=$minLen '{if ($NF > 99 && $4-$3 > M) print $1"\t"$6"\t"$2"\t"$7}' |sort |uniq > translation_hap1
cat assembly.mashmap.out |grep "$label2" |grep $g |awk -v M=$minLen '{if ($NF > 99 && $4-$3 > M) print $1"\t"$6"\t"$2"\t"$7}' |sort |uniq > translation_hap2

minLen=`echo $minLen |awk '{print $1*15}'`
cat translation_hap1|sort -k2,2|awk -v M=$minLen '{if ($3 > M) print $0}' |awk -v LAST="" -v S="" '{if (LAST != $2) { if (S > 0) print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG; SUM=0; MAX=0; C=0; } LAST=$2; S=$NF; SUM+=$3; if (MAX < $3) MAX=$3; C+=1; TIG=$1} END {print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG;}'  |awk '{print $1"\t"$4}' |sort -nk1,1 -s > chr_completeness_max_hap1

if [ -s translation_hap2 ]; then
   cat translation_hap2|sort -k2,2|awk -v M=$minLen '{if ($3 > M) print $0}' |awk -v LAST="" -v S="" '{if (LAST != $2) { if (S > 0) print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG; SUM=0; MAX=0; C=0; } LAST=$2; S=$NF; SUM+=$3; if (MAX < $3) MAX=$3; C+=1; TIG=$1} END {print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG;}'  |awk '{print $1"\t"$4}' |sort -nk1,1 -s > chr_completeness_max_hap2
fi

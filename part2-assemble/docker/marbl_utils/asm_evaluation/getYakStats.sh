#!/bin/bash

# given an assembly and some yak databases along with a ODB10 folder, this will produce some summary statistics for each haplotype (and combined in the case of a hic run)
# assumes there is a folder named *_odb10 in the current folder, for example primate_odb10 for human
# assumes there is a child.yak, mat.yak, and pat.yak files in the current folder
# will optinally also run asmgene if the ref.cdna.paf and reference cdna sequences are available, will be skipped if the file doesn't exist

yak=$(which yak 2>/dev/null)
if [ "x$yak" == "x" ]; then
   module load yak
   yak=$(which yak 2>/dev/null)
fi
if [ "x$yak" == "x" ]; then
   echo "Error: yak not found"
   exit
fi
minimap2=$(which minimap2 2>/dev/null)
if [ "x$minimap2" == "x" ]; then
   module load minimap2
   minimap2=$(which minimap2 2>/dev/null)
fi
if [ "x$minimap2" == "x" ]; then
   echo "Error: minimap2 not found"
   exit
fi
compleasm=$(which compleasm.py 2>/dev/null)
if [ "x$compleasm" == "x" ]; then
   compleasm="/data/korens/devel/compleasm_kit/compleasm.py"
fi
if [ "x$compleasm" == "x" ]; then
   echo "Error: compleasm not found"
   exit
fi
compleasm_db=`dirname $compleasm`
compleasm_db="$compleasm_db/mb_downloads"
db="primates_odb10"
if [[ "$#" -ge 1 ]]; then
   echo "Using custom db $1"
   db=$1
fi
if [ ! -d $compleasm_db/$db ]; then
   echo "Error, unable to find compleasm db $db in $compleasm_db, check download path/name and try again"
   exit
fi

if [ x$SLURM_CPUS_PER_TASK = x -o x$SLURM_CPUS_PER_TASK = xundefined -o x$SLURM_CPUS_PER_TASK = x0 ]; then
  cores=`grep -c ^processor /proc/cpuinfo`
else
  cores=$SLURM_CPUS_PER_TASK
fi

p=`pwd`


p=`dirname $p`
echo $p


if [ ! -e mat.yak ] && [ ! -L mat.yak ]; then
   echo "Error: missing mat.yak file"
   exit
fi
if [ ! -e pat.yak ] && [ ! -L pat.yak ]; then
   echo "Error: missing pat.yak file"
   exit
fi
if [ ! -e child.yak ] && [ ! -L child.yak ]; then
   echo "Error: missing child.yak file"
   exit
fi

for i in `ls assembly.haplotype[12].fasta`; do
   prefix=`echo $i |sed s/.fasta//g`

   if [ ! -e $prefix.yak.trioeval ] || [ $prefix.yak.trioeval -ot $i ]; then
      $yak qv -t $cores -l 100000 child.yak $i > $prefix.yak.qv
	  sleep 10
      $yak trioeval -t $cores mat.yak pat.yak $i > $prefix.yak.trioeval
	  sleep 10
   fi

   if [ ! -e $prefix.$db.full_table.tsv ] || [ $prefix.$db.full_table -ot $i ]; then
      $compleasm run -t$cores -l $db --library_path $compleasm_db -a $i -o $prefix
      mv $prefix/summary.txt $prefix.$db.summary.txt
      mv $prefix/$db/full_table.tsv $prefix.$db.full_table.tsv
   fi

   if [ -e ref.cdna.paf ]; then
      if [ ! -s $prefix.cdna.paf ]; then
         $minimap2 -cxsplice:hq -t$cores $i Homo_sapiens.GRCh38.cdna.all.fa > $prefix.cdna.paf
         paftools.js asmgene -i0.97 -a ref.cdna.paf $prefix.cdna.paf > $prefix.stats
      fi
   fi
done

if [ ! -e assembly.yak.trioeval ] || [ assembly.yak.trioeval -ot assembly.fasta ]; then
  $yak qv -t $cores -l 100000 child.yak assembly.fasta > assembly.yak.qv
  sleep 10
  $yak trioeval -t $cores mat.yak pat.yak assembly.fasta > assembly.yak.trioeval
  sleep 10
fi

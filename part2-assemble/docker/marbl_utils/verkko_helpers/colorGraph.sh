#!/bin/bash

color() {
   local input=$1
   local threshold=$2
   local mincount=100
   local prefix=`echo $input |sed s/.hapmers.count//g`

   cat $input | awk -v M=$mincount -v L=$threshold  'BEGIN \
   { \
     FS="[ \t]+"; OFS="\t"; \
     print "node\tmat\tpat\tmat:pat\tcolor"; \
   } \
   $1 != "Assembly" \
   { color = "#AAAAAA"; \
      if ($3+$4 > M) { \
         if      ($3 > ($3+$4)*L) { color = "#FF8888"; } \
         else if ($4 > ($3+$4)*L) { color = "#8888FF"; } \
         else                                        { color = "#FFFF00"; } \
      } \
      print $2, $3, $4, $3 ":" $4, color; \
   }'  > $prefix.colors.csv
}

if [ ! -e $MERQURY/trio/hap_blob.sh ]; then
   echo "Error: merqury not found, expected an env var named '\$MERQURY'"
   exit
fi

verkko=$(which verkko 2>/dev/null)
if [ "x$verkko" == "x" ]; then
   echo "Error: verkko not found"
   exit
fi
verkko=`dirname $verkko |awk '{print $1"/../lib/verkko/scripts/inject_coverage.py"}'`
if [ ! -e $verkko ]; then
   echo "Error: no verkko src found"
   exit
fi

asm=$1
mat=`cat $asm/verkko.yml |grep ruk_hap1 |awk '{print substr($NF, 2, length($NF)-2)}'`
pat=`cat $asm/verkko.yml |grep ruk_hap2 |awk '{print substr($NF, 2, length($NF)-2)}'`

if [ ! -e $asm/2-processGraph/unitig-unrolled-hifi-resolved.hifi-coverage.csv ]; then
   echo "Error: asm seems to be incomplete, missing 2-procesGraph output"
   exit 1
fi
if [ ! -e $asm/4-processONT/unitig-unrolled-ont-resolved.hifi-coverage.csv ]; then
   echo "Error: asm seems to be incomplete, missing 4-processONT output"
   exit 1
fi
if [[ ! -d $mat  ||  ! -d $pat ]]; then
   echo "Error: asm seems to be missing trio info, can't find one of $mat $pat"
fi
echo "Running trio on asm $asm and dbs are mat: $mat and pat: $pat"

cd $asm/2-processGraph
if [ ! -s unitig-unrolled-hifi-resolved.colors.csv ]; then
   cat unitig-unrolled-hifi-resolved.gfa  |awk '{if (match($1, "^S")) { print $1"\t"$2"\t*\tLN:i:"length($3)} else print $0}' > tmp.gfa
   $verkko --allow-absent unitig-unrolled-hifi-resolved.hifi-coverage.csv  tmp.gfa > unitig-unrolled-hifi-resolved.noseq.gfa
   cat unitig-unrolled-hifi-resolved.gfa  |awk '{if (match($1, "^S")) {print ">"$2; print $3}}'|fold -c > unitig-unrolled-hifi-resolved.fasta
   sh $MERQURY/trio/hap_blob.sh $mat $pat unitig-unrolled-hifi-resolved.fasta unitig-unrolled-hifi-resolved
   color unitig-unrolled-hifi-resolved.hapmers.count 0.90
fi

cd ../4-processONT
if [ ! -s unitig-unrolled-ont-resolved.colors.csv ]; then
   cat unitig-unrolled-ont-resolved.gfa  |awk '{if (match($1, "^S")) { print $1"\t"$2"\t*\tLN:i:"length($3)} else print $0}' > tmp.gfa
   $verkko --allow-absent unitig-unrolled-ont-resolved.hifi-coverage.csv  tmp.gfa > unitig-unrolled-ont-resolved.noseq.gfa
   cat unitig-unrolled-ont-resolved.gfa  |awk '{if (match($1, "^S")) {print ">"$2; print $3}}'|fold -c > unitig-unrolled-ont-resolved.fasta
   sh $MERQURY/trio/hap_blob.sh $mat $pat unitig-unrolled-ont-resolved.fasta unitig-unrolled-ont-resolved
   color unitig-unrolled-ont-resolved.hapmers.count 0.90
fi
cd ../..

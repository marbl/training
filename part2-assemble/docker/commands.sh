#!/bin/bash

# this script contains commands from assembly tutorial at AGBT-Ag
# they should work with the provided docker container with the exception of those which require large input files/input that have been omitted for efficiency
# note that this is not a complete assembly or data image/download, many of the outputs are just empty placeholders to provide runnable examples
# this is just give you an idea of what analysis can be run after an assembly is generated, presumably you will have your own data and assembly
# the data for this asembly came from the genomeark: https://www.genomeark.org/t2t-all/
# the initial asm was generated via
#verkko --screen human -d verkko_v2.0/  --hifi `pwd`/hifi/*.fastq.gz --nano `pwd`/ont/*.fast[aq].gz --hap-kmers maternal.k30.hapmer.meryl paternal.k30.hapmer.meryl trio

# overview of output files
cd data/verkko_v2.0
ls -lha assembly.*

# build translation between hifi and final nodes
python /marbl_utils/verkko_helpers/utig4_to_utig1.py . > utig4_2_utig1

# get t2t stats  
bash /marbl_utils/asm_evaluation/getT2T.sh assembly.fasta

# map to reference and generate names for chromosomes to plot
# mapping will take a bit of time, here we have pre-computed alignment files to save time
bash /marbl_utils/asm_evaluation/getChrNames.sh reference.fasta
ls *.csv
head translation_hap*

# remove human rDNA
bash /marbl_utils/verkko_helpers/removeRDNA.sh rDNA.fasta
ls *.gfa

# remove the candidate satellite sequence we found
# take the sequence of a single node with the candidate and extract it
grep utig4-354 *scfmap |awk '{print $2}' > tmp
seqtk subseq assembly.fasta tmp > repeat.fasta
# remove redundant sequence on the ends of repeat unit so it's only one copy
python /opt/conda/envs/myenv/lib/verkko/scripts/circularize_ctgs.py -f 0.75 repeat.fasta --min-ovl 1000 -o repeat.circular.fasta 
# make a temporary directory for our outputs and move/link the files we need to run the script to it
mkdir removeRepeat
cd removeRepeat
mv ../repeat.fasta ./
mv ../repeat.circular.fasta ./
ln -s ../assembly.colors.csv 
ln -s ../assembly.fasta
ln -s ../assembly.homopolymer-compressed.chr.csv 
ln -s ../assembly.homopolymer-compressed.noseq.gfa 
ln -s ../assembly.paths.tsv 
ln -s ../assembly.scfmap 
ln -s ../assembly.telomere.bed
ln -s ../6-rukki/
bash /marbl_utils/verkko_helpers/removeRDNA.sh repeat.circular.fasta 
cd ../..

# now run scaffolding using this repeat rather than the default of human rDNA
# this creates the necessary structure to run verkko re-using the current assembly and only mapping HiC data to it
# this command will take a while and requires all the data data that is not part of the tutorial
# it also requires the full assembly from v2.0 which again is not part of the tutorial
# the hifi should be in the folder hifi, ONT in folder ont, and hic in folder hic
# the data can be downloaded from https://www.genomeark.org/t2t-all/
# note it's always a good idea to run first with dry-run just to make sure the run will process what you expect, in this case it should start at 6-layout (e.g. after the linked files)
#mkdir verkko_hic_repeat
#cd verkko_hic_repeat
#ln -s  ../verkko_v2.0/1-buildGraph/
#ln -s  ../verkko_v2.0/2-processGraph/
#ln -s  ../verkko_v2.0/3-align
#ln -s  ../verkko_v2.0/3-alignTips/
#ln -s  ../verkko_v2.0/4-processONT/
#ln -s  ../verkko_v2.0/5-untip/
#cd ../
# verkko --screen human --rdna-scaff --snakeopts '--dry-run -U hicPhasing' --slurm -d verkko_hic_repeat/ --hifi  `pwd`/hifi/*.fastq.gz --nano `pwd`/ont/*.fast[aq].gz --hic1 hic/*R1*fastq.gz --hic2 hic/*R2*fastq.gz |less
# verkko --screen human --rdna-scaff --snakeopts '-U hicPhasing' --slurm -d verkko_hic_repeat/ --hifi  `pwd`/hifi/*.fastq.gz --nano `pwd`/ont/*.fast[aq].gz --hic1 hic/*R1*fastq.gz --hic2 hic/*R2*fastq.gz
# now ignore the HiC phasing and retain the original trio phasing
cd verkko_hic_repeat/8-hicPipeline
mv hicverkko.colors.tsv hicverkko.hiccolors.tsv
cp ../../verkko_v2.0/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.colors.csv ./hicverkko.colors.tsv
cp ../../verkko_v2.0/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.paths.gaf ./rukki.paths.gaf
cp ../../verkko_v2.0/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.paths.tsv ./rukki.paths.tsv
cd ../../
verkko --screen human --rdna-scaff --snakeopts '--touch   -U HiC_rdnascaff' --rdna-scaff-reference `pwd`/verkko_v2.0/removeRepeat/repeat.circular.fasta -d verkko_hic_repeat/  --hifi `pwd`/hifi/*.fastq.gz --nano `pwd`/ont/*.fast[aq].gz --hic1 hic/*R1*fastq.gz --hic2 hic/*R2*fastq.gz > /dev/null 2>&1
verkko --screen human --rdna-scaff --snakeopts '--dry-run -U HiC_rdnascaff' --rdna-scaff-reference `pwd`/verkko_v2.0/removeRepeat/repeat.circular.fasta -d verkko_hic_repeat/  --hifi `pwd`/hifi/*.fastq.gz --nano `pwd`/ont/*.fast[aq].gz --hic1 hic/*R1*fastq.gz --hic2 hic/*R2*fastq.gz |less
verkko --screen human --rdna-scaff --snakeopts '          -U HiC_rdnascaff' --rdna-scaff-reference `pwd`/verkko_v2.0/removeRepeat/repeat.circular.fasta -d verkko_hic_repeat/  --hifi `pwd`/hifi/*.fastq.gz --nano `pwd`/ont/*.fast[aq].gz --hic1 hic/*R1*fastq.gz --hic2 hic/*R2*fastq.gz
grep rdna verkko_hic_repeat/8-hicPipeline/rukki.paths.tsv

cd verkko_v2.0
# now moving on to resolving tangle examples
# these commands require the original 3-align/split folder which is not included in the download once again
# they are listed here for convience but do not run locally
#mkdir 8-manualResolution
cd 8-manualResolution
#touch empty.fasta
#GraphAligner -t 1 -g ../5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.gfa -f empty.fasta -a empty.gaf \
#	--diploid-heuristic 21 31 --diploid-heuristic-cache diploid.index
#	--seeds-mxm-cache-prefix $prefix \  --bandwidth 15 \
#	--seeds-mxm-length 30 \
#	--mem-index-no-wavelet-tree \
#	--seeds-mem-count 10000 && touch graph.index
#for i in `ls 3-align/split/*.fasta.gz`; do
#	id=`echo $i |sed s/.fasta.gz//g`
#	GraphAligner -t 24 -g unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.gfa -f $i -a $id.WORKING.gaf \
#		--seeds-mxm-windowsize 5000 --seeds-mxm-length 30 \
#		--seeds-mem-count 10000 \
#		--bandwidth 15 \
#		--multimap-score-fraction 0.99 \
#		--precise-clipping 0.85 \
#		--min-alignment-score 5000 \
#		--hpc-collapse-reads \
#		--discard-cigar \
#		--clip-ambiguous-ends 100 \
#		--overlap-incompatible-cutoff 0.15 \
#		--max-trace-count 5 \
#		--mem-index-no-wavelet-tree \
#	&& \
#	mv -f $id.WORKING.gaf $id.gaf
#done
#cat aligned*gaf > ont-aligned.gaf
#cd ..

# example tangle 1 utig4-751 with no clear coverage (0x) and composed of multiple utig1-nodes w/no clear coverage 
grep utig4-752 ../utig4_2_utig1
grep utig4-752 ../utig4_2_utig1 |awk '{print $2}'|tr '>' '\n' |tr '<' '\n' |sort |uniq > tmp
grep -w -f tmp ../2-*/*hifi-coverage.csv
grep utig4-752 ont-aligns.gaf |grep utig4-753

# example tangle 2 utig4-509
grep utig4-508 ont-aligns.gaf |grep utig4-505 |grep -c utig4-506                                                                                                                                                                                                                                                            
grep utig4-508 ont-aligns.gaf |grep utig4-505 |grep -c utig4-507                                                                                                                                                                                                                                                            

# example tangle 3 utig4-794
grep utig4-794 ont-aligns.gaf |grep utig4-797

# example tangle 4 utig4-1260
grep utig4-1260 ont-aligns.gaf |grep -c utig4-1262
grep utig4-1260 ont-aligns.gaf |grep -c utig4-717
grep utig4-1260 ont-aligns.gaf |grep -c utig4-622

# last example where we have a complex tangle
# do not run this command as above it is slow to run in the image
#GraphAligner -t 24 -g ../5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.gfa -f bH_ASM_new_with_unchanged_with_Y.fasta -a alignedASM.gaf \
#	--diploid-heuristic 21 31 --seeds-mxm-length 30 \
#	--seeds-mem-count 10000 \
#	--bandwidth 15 \
#	--multimap-score-fraction 0.99 \
#	--precise-clipping 0.95 \
#	--min-alignment-score 10000 \
#	--hpc-collapse-reads \
#	--discard-cigar \
#	--clip-ambiguous-ends 100 \
#	--overlap-incompatible-cutoff 0.15 \
#	--max-trace-count 5 \
#	--mem-index-no-wavelet-tree \
# the below is an example of using a custom sequence to patch an assembly, this could also be used with individual ONT reads
# once you find a candidate, you may want to shrink the reference target to a smaller region and re-map, which is what I did here
# make sure to adjust coordinates to non-HPC space (usually multiply by 1.5)
grep chr10_new alignedASM.gaf |sort -snk3,3
samtools faidx bH_ASM_new_with_unchanged_with_Y.fasta chr10_new:87396890-88888333 > subset.fasta
#GraphAligner -t 24 -g ../5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.gfa -f bH_ASM_new_with_unchanged_with_Y.fasta -a subset.gaf \
#  --diploid-heuristic 21 31 --seeds-mxm-length 30 \
#  --seeds-mem-count 10000 \
#  --bandwidth 15 \
#  --multimap-score-fraction 0.99 \
#  --precise-clipping 0.95 \
#  --min-alignment-score 10000 \
#  --hpc-collapse-reads \
#  --discard-cigar \
#  --clip-ambiguous-ends 100 \
#  --overlap-incompatible-cutoff 0.15 \
#  --max-trace-count 5 \
#  --mem-index-no-wavelet-tree \
cat subset.gaf
# remove lines that look spurious
/opt/conda/envs/myenv/lib/verkko/scripts/insert_aln_gaps.py ../assembly.homopolymer-compressed.gfa subset.gaf 1 50000 patch.nogap.gaf patch.gaf gapmanual y > patch.gfa

cp  ../6-layoutContigs/combined-alignments.gaf ./
cat patch.gaf >> combined-alignments.gaf

cp  ../6-layoutContigs/combined-edges.gfa ./
cat patch.gfa | grep '^L' |grep gap >> combined-edges.gfa

ln -s ../6-layoutContigs/combined-nodemap.txt
cp ../6-layoutContigs/nodelens.txt ./
cat  patch.gfa | grep gap | awk 'BEGIN { FS="[ \t]+"; OFS="\t"; } ($1 == "S") && ($3 != "*") { print $2, length($3); }' >> nodelens.txt
cp ../7-consensus/ont_subset.fasta.gz ./
cat subset.fasta | gzip -c >> ont_subset.fasta.gz
seqtk gc ont_subset.fasta.gz |awk '{print $1}'|sort |uniq > ont_subset.id

cp ../../verkko_hic_repeat/8-hicPipeline/rukki.paths.gaf ./manual.paths.gaf
cat patch.gaf
# manual edit the file to add the paths we got above from our traversals and also the gapfilled one above

# this will confirm that the layout is valid, once done, you can make a new assembly folder with the links to the existing 1-5 folders, a new 6-layoutContigs file
# the command will likely not fit in docker so don't expect it to fully complete 
/opt/conda/envs/myenv/lib/verkko/scripts/get_layout_from_mbg.py combined-nodemap.txt combined-edges.gfa combined-alignments.gaf manual.paths.gaf nodelens.txt unitig-popped.layout unitig-popped.layout.scfmap
cd ../../

mkdir verkko_final_asm
cd verkko_final_asm
ln -s  ../verkko_v2.0/1-buildGraph/
ln -s  ../verkko_v2.0/2-processGraph/
ln -s  ../verkko_v2.0/3-align
ln -s  ../verkko_v2.0/3-alignTips/
ln -s  ../verkko_v2.0/4-processONT/
ln -s  ../verkko_v2.0/5-untip/
mkdir 6-layoutContigs
cd 6-layoutContigs
ln -s ../../verkko_v2.0/8-manualResolution/combined-nodemap.txt
ln -s ../../verkko_v2.0/8-manualResolution/combined-edges.gfa
ln -s ../../verkko_v2.0/8-manualResolution/combined-alignments.gaf
ln -s ../../verkko_v2.0/8-manualResolution/nodelens.txt
cd ..
mkdir 6-rukki
cd 6-rukki
ln -s ../../verkko_v2.0/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.noseq.gfa
ln -s ../../verkko_v2.0/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.colors.csv
ln -s ../../verkko_v2.0/8-manualResolution/manual.paths.gaf unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.paths.gaf
ln -s ../../verkko_v2.0/8-manualResolution/manual.paths.gaf unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.paths.tsv
cd ../
mkdir 7-consensus
cd 7-consensus
ln -s ../../verkko_v2.0/8-manualResolution/ont_subset.id
ln -s ../../verkko_v2.0/8-manualResolution/ont_subset.fasta.gz
cd ../
cd ../
verkko --screen human --snakeopts '--touch'   -d verkko_final_asm/  --hifi `pwd`/hifi/*.fastq.gz --nano `pwd`/ont/*.fast[aq].gz --hap-kmers maternal.k30.hapmer.meryl paternal.k30.hapmer.meryl trio > /dev/null 2>&1
verkko --screen human --snakeopts '--dry-run' -d verkko_final_asm/  --hifi `pwd`/hifi/*.fastq.gz --nano `pwd`/ont/*.fast[aq].gz --hap-kmers maternal.k30.hapmer.meryl paternal.k30.hapmer.meryl trio |less
#verkko --screen human                        -d verkko_final_asm/  --hifi `pwd`/hifi/*.fastq.gz --nano `pwd`/ont/*.fast[aq].gz --hap-kmers maternal.k30.hapmer.meryl paternal.k30.hapmer.meryl trio

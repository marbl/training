# Post assembly curation

Please refer to [Part 3 - Slides](https://docs.google.com/presentation/d/1lFy3Iy-po-GrKbqtx_nSYcT4R4f51CVxzM_q69Ag2VY/edit?usp=drive_link).

## Requirements
* [T2T-Polish](https://github.com/arangrhie/T2T-Polish): latest github tip
* [merqury](https://github.com/marbl/merqury): latest github tip
* [meryl](https://github.com/marbl/meryl): v1.4.1 release
* [seqrequester](https://github.com/marbl/seqrequester): latest github tip
* [seqtk](https://github.com/lh3/seqtk)
* [Winnowmap2](https://github.com/marbl/winnowmap) 2.03 release
* [minimap2 and k8](https://github.com/lh3/minimap2) for sam2paf.js. See [misc](https://github.com/lh3/minimap2/tree/master/misc)
* [asset](https://github.com/dfguan/asset/)
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [ucsc utils](https://hgdownload.soe.ucsc.edu/downloads.html#source_downloads): wigToBigWig
* java
* R

This tutorial assumes T2T-Polish, seqrequester, Winnowmap, asset are installed under `$tools`. The rest are assumed to be set in your environment $PATH (`~/.bash_profile`) or loaded via modules.

## 1. Clean up

1. Identify MT, chloroplast, or any other contamination

  - Mash screen and Mash distance: https://github.com/VGP/vgp-assembly/tree/master/pipeline/mash
  - NCBI FCS: https://github.com/ncbi/fcs

2. Remove excessive small copies of rDNA, small repetitive sequences, ...

3. Add ribotin rDNA (if available) and MT sequences

4. Name chromosomes (eg. chr1_1, chr1_2, ...)

## 2. Coverage based analysis

### Telomere, gaps, microsatellite annotation
```shell
mkdir -p pattern && cd pattern
asm=bH_v2
ver=v2    # Tip: I prefer the polished ones be v1.0 or v2.0. Set accordingly
ln -s ../assemblies/$asm.fasta
$tools/T2T-Polish/coverage/init.sh $asm.fa $ver
cd ../
```

Output files:
```
  $ver.bed          # size
  $ver.exclude.beds # N-gaps
  $ver.telo.bed     # telomere
  # Microsatellite repeats (pattern)
  $ver/$ver.microsatellite.AT.128.bw
  $ver/$ver.microsatellite.GA.128.bw
  $ver/$ver.microsatellite.GC.128.bw
  $ver/$ver.microsatellite.TC.128.bw
  # ... (and .bed for the above bw files)
```
Load them on IGV and investigate satellite repeat enriched region. Log along scaffolds that should be re-oriented to match the p and q arm definitions and/or match between haplotypes. Re-orienting can happen after read-mapping and polishing.

### Mapping reads back to the assembly
Align all reads to 
* both haplotypes
* each haplotype

As a result, there will be 3 x 3 alignments:<br>
(HiFi, ONT, Illumina) x (Both haps, hap1, hap2)
```shell
mkdir -p mapping && cd mapping
ls /full/path/to/HiFi/*.fq.gz > HiFi.fofn
ls /full/path/to/ONT/*.fq.gz > ONT.fofn
ls /full/path/to/Illumina/*read-pair1.fq.gz > illumina1.fofn
ls /full/path/to/Illumina/*read-pair2.fq.gz > illumina2.fofn
paste illumina1.fofn illumina2.fofn > illumina.fofn
rm illumina1.fofn illumina2.fofn
```
For `slurm` environment, use the submitter wrapper. Otherwise, run each scripts described as in the [README](https://github.com/arangrhie/T2T-Polish/tree/master/winnowmap).
Understand winnowmap submitter script
```shell
$tools/T2T-Polish/winnowmap/_submit.sh
Usage: ./_submit.sh ref.fasta prefix map [wm_opt]
  ref.fasta  reference to align
  prefix     output prefix
  map        mapping mode. map-pb for HiFi, map-ont for ONT, map-pb-clr for CLR.
  wm_opt     winnowmap optional arguments. i.e. -y for adding methylation tags.
  Required: input.fofn
```

For example, for hifi reads, submit for the 3 assemblies:
```shell
pf=HiFi
asm=bH_v2 # adjust
asm_fa=../assemblies/$asm
# For convenience, we assume ../assemblies/ has all 3 versions as $asm_fa.$hap.fa

for hap in dip hap1 hap2
do
  mkdir -p ${pf}_to_$hap && cd ${pf}_to_$hap
  ln -s ../HiFi.fofn input.fofn
  $tools/T2T-Polish/winnowmap/_submit.sh $asm_fa.$hap.fa ${pf}_to_$hap map-pb
  cd ../
done
```
Submit the same script for ONT reads; using `ONT.fofn` as `input.fofn` and replacing `map-pb` to `map-ont`.

Illumina reads can be mapped with `bwa`, as in [this script](https://github.com/arangrhie/T2T-Polish/tree/master/bwa).

## 3. k-mer based errors and haplotype switches
Let's get back to the upper directory from mapping and build kmer dbs for merqury. Re-use the `input.fofn` files in mapping. In this example, we will generate k=31 mers for evaluation.

```shell
cd ../
mkdir -p meryl && cd meryl
for pf in illumina HiFi ONT
do
  mkdir $pf && cd $pf
  ln -s ../../mapping/$pf.fofn
  $MERQURY/_submit_build.sh 31 $pf.fofn $pf # with slurm
  cd ../
done
```
Once the meryl dbs are built, they will appear as `$pf.k31.meryl`. Sym-link them all under meryl for convenience.
```shell
cd ../
for pf in illumina HiFi
do
  ln -s $pf/$pf.k31.meryl $pf.meryl
done
```

### For trios
Apply the same for parental datasets; to obtain e.g. `pat.k31.meryl` and `mat.k31.meryl` if available. Sym-link as the rest of the meryl dbs.

Generate [hapmers](https://github.com/marbl/merqury/wiki/1.-Prepare-meryl-dbs#3-build-hap-mer-dbs-for-trios):
```shell
cd ../
mkdir -p hapmers && cd hapmers
$MERQURY/hapmers.sh ../mat.meryl ../pat.meryl ../illumina.meryl
cd ../
```
### Hybrid db
In addition to the raw reads kmer sets, we need the hybrid kmer db - which is a merged version of illumina and HiFi kmers.
See [T2T-Polish](https://github.com/arangrhie/T2T-Polish/tree/master/merqury) for more details.

```shell
meryl union-sum [ greater-than 1 illumina.k31.meryl ] [ greater-than 1 HiFi.k31.meryl ] output hybrid.meryl
```

### Merqury
Now we are ready to generate kmer profiles.
```shell
cd ../
mkdir -p merqury && cd merqury
```

Run Merqury on each platform: `illumina`, `HiFi`, `hybrid`.
```shell
ln -s ../meryl/illumina.meryl
ln -s ../meryl/hapmers/mat.hapmer.meryl mat.meryl
ln -s ../meryl/hapmers/pat.hapmer.meryl pat.meryl

# For trio based evaluation
mkdir -p illumina && cd illumina
$MERQURY/_submit_merqury.sh ../illumina.meryl ../mat.meryl ../pat.meryl ../../assemblies/$asm.hap1.fa ../../assemblies/$asm.hap2.fa illumina
cd ../

# HiFi and Hybrid - for QV estimates
for pf in HiFi hybrid # illumina if no hapmers are available
do
  ln -s ../meryl/$pf.meryl
  mkdir -p $pf && cd $pf
  $MERQURY/_submit_merqury.sh ../$pf.meryl ../../assemblies/$asm.hap1.fa ../../assemblies/$asm.hap2.fa $pf
  cd ../
done
```
Once Merqury completes, collect `*_only.bed` and `*_only.wig`.
```shell
cd hybrid
cat *_only.bed | bedtools merge -i - > $asm.error.bed
cd ../
```

Additionally, collect hapmer tracks for trios:

```shell
cd illumina
# Merge hapmer tracks
cat *.mat.wig > $asm.mat.wig
cat *.pat.wig > $asm.pat.wig

# Convert to bw files (optional, to reduce size)
wigToBigWig $asm.mat.wig ../../assemblies/$asm.fa.fai $asm.mat.bw
wigToBigWig $asm.pat.wig ../../assemblies/$asm.fa.fai $asm.pat.bw
cd ../
```

Link the $asm.error.bed under pattern.
```shell
cd ../pattern/
ln -s ../../merqury/hybrid/$asm.error.bed $asm.error.bed
cd ../
```

## Issues track

Winnowmap submitter generates filtered bam and paf files at the end. The most informative issues comes from read_to_dip alignments.

This step requires a lot of files to be in place. Make sure they exist in the proper location. Again, don't forget to link `$asm.error.bed` from Merqury hybrid.

```shell
$tools/T2T-Polish/coverage/issues.sh
Usage: ./issues.sh in.paf name ver platform pattern
  in.paf   input paf file
  name     name to appear on the .wig files
  ver      assembly version
  platform HiFi, ONT or Hybrid
  pattern  path to pattern folder

Required: pattern/ver/, pattern/ver.bed, pattern/ver.error.bed, pattern/ver.exclude.bed, pattern/ver.telo.bed
```
The above script generates a lot of `.wig` files and `.bed` files, as well as useful coverage related statistics.
Refer to the [coverage](https://github.com/arangrhie/T2T-Polish/tree/master/coverage) description for more details. `ver` should match as it appears in the `pattern`.

```shell
for pf in HiFi ONT
do
  out=${pf}_to_dip
  mkdir -p $out && cd $out
  paf=`ls ../../mapping/$out/$out.*pri.paf`
  ln -s $paf
  paf=`basename $paf`
  $tools/T2T-Polish/coverage/issues.sh $paf $out $ver $pf ../../pattern
  cd ../
done
```

### Evaluate on IGV
Load the `*.issues.bed` and `*.cov.wig` on IGV, along with the pattern `*.bw` files, and enjoy(?) reading the error profiles!

## Next steps
Use the alignments to call variants, inspect, and create correction `.vcf` files.
Alternatively, the alignments could get further corrected using [SeqPhase](https://github.com/mobinasri/secphase) and polished with [DeepPolisher](https://github.com/google/deeppolisher). Note that DeepPolisher currently only works for HiFi reads.
### Software versions
```
bwa 0.7.17
deeptools 3.3.0
picard 2.20.3 
seqkit 0.12.0
sambamba 0.7.0
java 1.6.0
samtools 1.2
htslib 1.2.1
bedtools 2.29.0
macs2 2.1.1.20160309
finder 2
Jaguar 1.7.5
```

### Pseudocode of pipelines

#### Alignment of short-reads

```
# index genome
bwa index reference_genome.fa

# align reads, convert sam to bam file and sort by coordinates
bwa mem -t 24 bwa_genome_index file_R1.fastq file_R2.fastq | sambamba view -S -h -f bam -t 24 /dev/stdin | sambamba sort -t 24 --tmpdir=/path/ /dev/stdin -o /outpath/file.sorted.bam

# mark duplicated reads
java -jar -Xmx10G MarkDuplicates.jar I=file.sorted.bam O=file.sorted.dups_marked.bam M=dups AS=true VALIDATION_STRINGENCY=LENIENT QUIET=true
```
#### ChIP-seq data analysis

```
## Peak calling with findER2 (in-house)
# Single-end libs 
java -jar -Xmx25G /home/mbilenky/ew/Solexa_Java/jarsToDeploy/finder2.jar inputBam:$inp signalBam:$sig outDir:$out acgtDir:/projects/epigenomics2/users/mbilenky/resources/hg38/ACGT SE

# Paired-end libs
$java -jar -Xmx25G /home/mbilenky/ew/Solexa_Java/jarsToDeploy/finder2.jar inputBam:$inp signalBam:$sig outDir:$out acgtDir:/projects/epigenomics2/users/mbilenky/resources/hg38/ACGT

## RPKM normalization of bw files
# Single-end libs 
for f in files*bam; 
do echo $f;
bamCoverage -b $f -o ../bigwig_RPKM_norm_50bp/$f.bw -of bigwig -bs 50 --effectiveGenomeSize 2913022398 --normalizeUsing RPKM --extendReads 175 --ignoreDuplicates -p max/2;
done

# Paired-end libs
for f in files*bam; 
do echo $f;
bamCoverage -b $f -o ../bigwig_RPKM_norm_50bp/$f.bw -of bigwig -bs 50 --effectiveGenomeSize 2913022398 --normalizeUsing RPKM --extendReads --ignoreDuplicates -p max/2;
done

# HOMER motif finding
findMotifsGenome.pl regionOfInterest.bed hg38 ./ size 200 -len 8 -p 16  

## Generate profile plot my merging samples
cd /path/bw/
m=ENS_RPKM_more0.5.v2.bed # protein conding genes with >=0.05 RPKM
l=ENS_RPKM_less0.5.v2.bed # protein conding genes with <0.05 RPKM

out=/path/

for bw in *CLL*bw; do echo $bw;

computeMatrix reference-point --referencePoint center \
                              --regionBodyLength 5000 \
                              --upstream 500 \
                              --downstream 500 \
                              -S  $bw \
                              -R $m $l \
                              -p 40 \
                              --smartLabels \
                              --skipZeros -o $out/deeptools/$bw.matrix.gz
done 

cd $out/deeptools/
computeMatrixOperations rbind -m *matrix.gz -o DNAme_moreLessgenes_merged_20CLL.gz
rm *.matrix.gz

plotProfile -m DNAme_moreLessgenes_merged_23CLL.gz \
     -out $out/deeptools/DNAme_moreLessgenes_merged_23CLL.pdf \
      --samplesLabel label_1 \
      --plotFileFormat 'pdf'

## Making tracks for UCSC genome browser
                                
# Example for H3K27ac                               
cd /gsc/edcc_web/CLL_rislam/H3K27ac/hg38/
                                
mark=H3K27ac
for i in 1;
do
echo "track $mark"
echo "compositeTrack on"
echo "shortLabel $mark"
echo "longLabel $mark"
echo "type bigWig"
echo ""

for f in *NBC*bw *GCBC*bw *MBC*bw *PBC*bw *uCLL*bw *mCLL*bw *_CLL_*;  
do 
echo $f | awk '{gsub(".bw", "") ; print $0}' | tr -s '_' '\t' | awk '{print "        " " track", $1"_"$4"_"$5"_"$6"_"$7$8}'
echo $f | awk '{gsub(".bw", "") ; print $0}' | tr -s '_' '\t' | awk '{print "        " " shortLabel", $5"_"$6"_"$7$8}'
echo $f | awk '{gsub(".bw", "") ; print $0}' | tr -s '_' '\t' | awk '{print "        " " longLabel", $1"_"$4"_"$5"_"$6"_"$7$8}'
echo "        " "parent $mark on"
echo "        " "type bigWig"
echo "        " "visibility full"
echo "        " "maxHeightPixels 70:70:32"
echo "        " "configurable on"
echo "        " "autoScale on"
echo "        " "alwaysZero on"
echo "        " "priority 0.1"
echo "        " "bigDataUrl $f"
echo "        " "color 0,0,255"
echo "        " ""
done 
done >trackDb.txt
```

#### RNA-seq data analysis

```
# repositioning by Jaguar. Every different read length has a different Jaguar reference! 
jaguar_ref=/path/ref.fa

# NOTE: BAM FOR JAGUAR MUST BE NAME SORTED!
bwa mem -M -P -t 24 $jaguar_ref <FQ1> <FQ2>  | sambamba_v0.5.5 view -S -h -f bam -t 16 /dev/stdin | sambamba_v0.5.5 sort -t 16 -n /dev/stdin -o <BAM>

# Run repositioning by Jaguar. RunJR.sh is an in-house script
# Example of command with hg19v69
RunJR.sh $out/$lib".sortedByName.bam" $out/j hg19_ens69

# RNAseq master (in-house) 
# RNAseqMaster.sh is the master script to run number of tools and generate an extensive QC report and RPKM matrix.
RNAseqMaster.sh <1: bam file (long path)> <2: name> <3: folder with all output (will be this PATH/name)> <4: species (hg19v66/...)> <5: strand specific(S)/regular(R)> <6: quality threshold> <7: running mask COVERAGE,RPKM,LEAKAGE,PROFILE,REPORT (1==run, 0==don't run))> <8: path to resource folder, e.g. /project/epigenomics/resources/> [<9: java>] [<10: java max heap space>
 
OUTPUT: coverage files and coverage distributions; note that for strand specific RNA-seq coverages are calculated for proper strand and then cat together
```

#### WGBS data analysis

```
## Bisulfite Conversion Rate Computation Using Novomethyl
# Novomethyl requires two preliminary steps: split of the bam file into CT and GA strands, samtools mpileup (steps 1 and 2 below)
samtools view -h <SortedDupedBAMfile> | grep -e "^@\|ZB:Z:CT" | $samtools view -ubS - | $samtools sort - <FileName>.<species>.CT
samtools view -h <SortedDupedBAMfile> | grep -e "^@\|ZB:Z:GA" | $samtools view -ubS - | $samtools sort - <FileName>.<species>.GA

## Running samtools mpileup
samtools mpileup -BC 0 -q 30 -f <PathToReferenceGenomeFastaFile> <FileName>.<species>.CT.bam <FileName>.<species>.GA.bam | gzip -c > <FileName>.<species>.mpileup.txt.gz

## Running Novomethyl (for human and lambda spike in control)
gunzip -c <FileName>.<species>.mpileup.txt.gz | <PathToNovomethyl> -o Consensus -% 2> <FileName>.<species>.methylation.log |gzip -c > <FileName>.<species>.Cmethyl.cons.bed.gz

## The last lines of the methylation log file shows the estimated bisulfite conversion rate.
tail -2 output_filename.methylation.log

## run Novo5mc 
java -jar -Xmx10G Novo5mC.jar -bam file.bam -out /path/  -genome $human_ref -q5 -F 1540 -minCoverage 3 -name output_filename -samtools samtools -regions 1 > /desired/path/to/your/log/file.log

OUTPUT: This generates a file bearing the information of number of reads were conver or unconverted follwoing bisulfite conversion.

## calculate avg methylation of genomic features

bed=H3K27ac_CLL_enriched_over_others_4col.bed
map=/path/map_5col_CpG_files/ #wgbs

for cpg in $map/*.CpG.bed; 
do
    echo $cpg;
    outfile=$(basename $cpg .bed);
    echo $outfile;
    echo $bed;
    bedmap --echo --mean $bed $cpg | awk '{print $4}' | awk -F"|" '{print $2}' >$outfile.dname;
done 

# add rowname
less $bed | awk '{print $1"_"$2"_"$3}' >A.dname

# make table
((echo *.dname |tr ' ' '\t') && (paste *.dname)) >H3K27ac_CLL_enriched_over_others_dname.tsv

rm *.dname
```
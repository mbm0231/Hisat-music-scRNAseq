# Mapping next-generation sequencing reads from Equus caballus intestinal organoids treated with Rotavirus using Hisat2 to Equcab3 index. 
## hisat2-build builds a HISAT2 index from a set of DNA sequences from Equcab3.
```
#!/bin/bash

#SBATCH --time 24:00:00
#SBATCH --job-name=Hisat2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=jumbo
#SBATCH --mem=512GB
#SBATCH --mail-type ALL
#SBATCH --mail-user mbmo231@uky.edu,farman@uky.edu
#SBATCH -A cea_farman_s24cs485g
#SBATCH -o /project/farman_s24cs485g/mbmo231/alignments/%x_%J.out
#Hisat2 image location

# establishing a variable with the singularity contaiiner that you need for Hisat2
Hisat2image=/share/singularity/images/ccs/conda/amd-conda2-centos8.sinf
# Same thing as above, but the whole command for Samtools
SAMTOOLS='singularity run --app samtools113 /share/singularity/images/ccs/ngstools/samtools-1.13+matplotlib-bcftoools-1.13.sinf samtool$
# establish a variable that we use as the index name
indexname='GCF_002863925.1_EquCab3.0_genomic'
# establish where the actual fasta file is
fasta=${indexname}'.fna.gz'
# establish the directory with the index inside of it
indexdir=/project/farman_s24cs485g/mbmo231/Horse_ref/${indexname}_index
# run the Hisat indexing on the fasta file that we talked about above
#singularity run --app hisat2221 $Hisat2image hisat2-build --large-index -p 16 GRCh38p14_index/GRCh38p14.fasta GRCh38p14_index/GRCh38p
singularity run --app hisat2221 $Hisat2image hisat2-build --large-index -p 16 /project/farman_s24cs485g/mbmo231/Horse_ref/GCF_002863925.1_EquCab3.0_genomic.fna.gz 

```
## Mapping next-generation sequencing reads using Hisat2 to Equcab3 index 

```
#!/bin/bash

#SBATCH --time 24:00:00
#SBATCH --job-name=Hisat2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=jumbo
#SBATCH --mem=512GB
#SBATCH --mail-type ALL
#SBATCH --mail-user mbmo231@uky.edu,farman@uky.edu
#SBATCH -A cea_farman_s24cs485g
#SBATCH -o /project/farman_s24cs485g/mbmo231/alignments/%x_%J.out
#Hisat2 image location

# establishing a variable with the singularity contaiiner that you need for Hisat2
Hisat2image=/share/singularity/images/ccs/conda/amd-conda2-centos8.sinf
# Same thing as above, but the whole command for Samtools
SAMTOOLS='singularity run --app samtools113 /share/singularity/images/ccs/ngstools/samtools-1.13+matplotlib-bcftoools-1.13.sinf samtools'
# establish a variable that we use as the index name
indexname='GCF_002863925.1_EquCab3.0_genomic'
# establish where the actual fasta file is
fasta=${indexname}'.fna.gz'
# establish the directory with the index inside of it
indexdir=/project/farman_s24cs485g/mbmo231/Horse_ref/${indexname}_index
# run the Hisat indexing on the fasta file that we talked about above
#singularity run --app hisat2221 $Hisat2image hisat2-build --large-index -p 16 GRCh38p14_index/GRCh38p14.fasta GRCh38p14_index/GRCh38p14
# for some reason establish f as the directory where your raw reads are
f=/project/farman_s24cs485g/mbmo231/RNAseqRota
# make that directory... [ ! -d $f ] && mkdir $f ... that would be a better way of doing that 
# it adds logic saying "if a directory $f does NOT exist (the ! infront of everything) > then mkdir $f
# alternatively, mkdir -p kinda does the same thing
#mkdir $f
# $dir variable?!? 
dir=/project/farman_s24cs485g/mbmo231
# copy the raw reads into /project/farman_s24cs485g/mbmo231/rawReads - R1
cp $dir/1_S27_R1_001.fastq.gz   $f/
# copy the raw reads into /project/farman_s24cs485g/mbmo231/rawReads - R2
cp $dir/1_S27_R2_001.fastq.gz   $f/

# do the actual aligning and then send everything through samtools for sorting (and convertion from SAM to BAM format)
singularity run --app hisat2221 $Hisat2image hisat2 \
-p 16 -x $indexdir/$indexname -1 $f/1_S27_R1_001.fastq.gz -2 $f/1_S27_R2_001.fastq.gz     \
--summary-file alignments/JI_1_S27_.txt --dta-cufflink   \
| $SAMTOOLS sort - -@ 16 -O bam -o alignments/HS_JI_accepted_hits_1_S27_.bam
```

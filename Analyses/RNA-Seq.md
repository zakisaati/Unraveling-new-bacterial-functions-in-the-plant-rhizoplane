# RNA-Seq analyses

Here we show the codes used for analyzing RNA-Seq data and search for differentially expressed genes. 

## Mapping

The first step needed is to align each sample with quality filtered reads against the bacterial genome. We used [Bowtie2](https://github.com/BenLangmead/bowtie2) for this.

Before performing the alignment, it is needed to create a reference database with the genome:

~~~
$ bowtie2-build Genome_contigs.fa CDVBN10_build/CDVBN10.build
~~~

Then, the alignment:

~~~
$ bowtie2 -p 36 -x CDVBN10_build/CDVBN10.build -q sample_X_qc_reads.fq -S sample_X.sam
~~~

## Formating SAM files

The following steps will need to use sorted and indexed BAM files instead of those in SAM format. Thus, we used [SamTools](https://github.com/samtools/samtools) to modify the files into the desired format.

~~~
$ for file in *.sam
do
replicate=${file%%.sam}
echo "samtools view -Sb ${file} > ${replicate}.bam"
echo "rm -rf ${file}"
done > sam_to_bam.sh
~~~

~~~
$ bash sam_to_bam.sh
~~~

~~~
$ samtools sort replicate.bam > replicate_sorted.bam
~~~
~~~
$ samtools index ::: replicate_sorted.bam
~~~

## Counting transcripts per CDS

We will use DESeq2 later for differential expression analysis. This program requires to provide a count of transcripts for each CDS. 
To create those counts, we ran [FADU](https://github.com/IGS/FADU), a quantification tool designed specifically for prokaryotic RNA-Seq analyses.

The input consisted on the General Feature Formatted (GFF) genome from RAST annotation and the files obtained above.

~~~
$ mkdir replicate_X_FADU_count/
~~~

~~~
$ julia fadu.jl -g genome.gff -b "replicate_sorted.bam" -o "replicate_X_FADU_count" -s "yes" -f "CDS"
~~~

There will be an ouput file for each sample or replicate. This file would have several columns. Let's see as an example some lines of this output (copied from [here](https://github.com/IGS/FADU))

~~~
featureID   uniq_len    num_alignments  counts  tpm
cds0    1017    4.5 2.81    1744.38
cds1    1194    3.5 2.48    1310.50
cds10   1161    0.0 0.00    0.00
cds100  591 0.0 0.00    0.00
cds1000 741 8.0 7.46    6358.81
cds1001 850 5.0 4.45    3310.38
cds1002 829 1.0 0.48    363.86
cds1003 1167    0.0 0.00    0.00
cds1004 1164    0.0 0.00    0.00
cds1005 816 0.0 0.00    0.00
~~~

We manually retained just the first (featureID) and the fourth (counts) columns. Then, we rounded the counts and create a folder with each of the files with rounded counts.

## Differential Expression Analysis

We used [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), in R environment, to find up- or down-regulated genes in the studied conditions.

The input consisted on the previously obtained rounded count files. It's important to say that these counts should not be normalized.

~~~
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design= ~ Condition)
~~~

~~~
dds <- DESeq(ddsHTSeq)
~~~
~~~
$ res <- results (dds)
~~~
~~~
$ write.csv(as.data.frame(res), file="/path/to/output_DESeq2.csv")
~~~

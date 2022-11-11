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



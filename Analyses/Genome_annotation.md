# Genome annotation

Here we show the codes used for the genome and its encoded proteome annotations. 

NOTE: we just include those that requires command-line usage, not those annotations made through web interface programs, such as RAST or EggNOG-Mapper.

## CAZys

The search for [CAZys](http://www.cazy.org/) (Carbohydrate Active EnzYmes) was done through the [run_dbcan.py](https://github.com/linnabrown/run_dbcan) command, which is a standalone tool for [dbCAN2](https://bcb.unl.edu/dbCAN2/).

The input for this analysis consisted on the [encoded proteome](Source_data/CDVBN10.faa) (fasta amino acid file).

~~~
$ run_dbcan.py CDVBN10.faa protein --out_dir CDVBN10_dbcan2 --dia_cpu 36 --hmm_cpu 36 --hotpep_cpu 36 --tf_cpu 36 --tf_cpu 36 --db_dir /path/to/db
~~~

Then, as recoomended by authors for dbCAN2, we just retained those CAZYs that have been annotated through 2 of the 3 search methods. 

~~~
$ awk '$5 ~ /[23]/ { print $0 }' CDVBN10_dbcan2/overview.txt > CAZYs_2_or_3_tools.txt
~~~

## Biosynthetic Gene Clusters (BGCs)

We ran [antiSMASH](https://docs.antismash.secondarymetabolites.org/) (v5.1) to find Biosynthetic Gene Clusters (BGCs) related with the biosynthesis of secondary metabolites.

The RAST annotated genome in *Gene BanK* format was the [input](Source_data/CDVBN10.gbk) for this annotation.

~~~
antismash --fullhmmer --cb-general --cb-subclusters --cb-knownclusters --asf --pfam2go --genefinding-tool none --output-dir antismash_output CDVBN10_genome.gbk -c 36
~~~

## Signal Peptides


To look for proteins with a signal peptide we used the **SignalP** tool. We downloaded it from: https://services.healthtech.dtu.dk/cgi-bin/sw_request

Again, the input for this analysis consisted on the encoded proteome (fasta amino acid file). The command used was:

~~~
$ signalp -fasta CDVBN10.faa -org gram- -format short -prefix signalp_output.faa
~~~


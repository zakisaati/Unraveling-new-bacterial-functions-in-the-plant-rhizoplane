# Metagenomic screening for plant-associated genes and yafL

This section describes the methods used to screen metagenomes for homologs of plant-associated proteins and the phage gene yafL across multiple publicly available datasets.

## Preparation of training files and databases

Before launching the analysis, a copy of the training directory used by FragGeneScan was made locally to avoid path-related errors:

~~~
$ cp -r ~/miniconda3/envs/sra-toolkit/bin/train/ ./
~~~

Protein sequences for plant-associated genes and yafL were then compiled and used to build DIAMOND databases for downstream similarity searches:

~~~
$ diamond makedb --in up_L2fc_p0.05.faa -d plant_associated_db
$ diamond makedb --in yafL.faa -d yafL_db
$ pigz *.faa
~~~

## Workflow execution across multiple metagenomes
A custom Bash script was created to automate the analysis across multiple SRA/ENA accessions. The script was written into a file named workflow_metagenomic_screening.sh

~~~
$ nano workflow_metagenomic_screening.sh
~~~

The following code was then pasted and executed:

~~~
#!/bin/bash

# Output header
echo -e "SampleID\tPlantHits_85\tYafLHits_85\tPlantHits_75\tYafLHits_75\tTotalReads\tPlantHitsPerMillion_85\tYafLHitsPerMillion_85\tPlantHitsPerMillion_75\tYafLHitsPerMillion_75" > combined_results_normalized.txt

# List of samples to process. Insert as many SRA accession IDs (eg. SRR30765289 SRR30765293 SRR30765313 etc) as needed.
for id in SRRXXXXXXX SRRYYYYYYY
do
  echo "Proccessing sample: $id"

  # Download and extract FASTQ
  prefetch "$id"
  fasterq-dump "$id" --split-files --min-read-len 75 -p -e 6 --skip-technical

  # Quality filtering
  fastp -i "${id}_1.fastq" -I "${id}_2.fastq" \
        -o "${id}_1.clean.fastq" -O "${id}_2.clean.fastq" \
        --detect_adapter_for_pe --length_required 75 \
        --html "${id}_fastp_report.html" --thread 6

  rm -f "${id}_1.fastq" "${id}_2.fastq"

  # Truncate to 75 bp
  usearch -fastx_truncate "${id}_1.clean.fastq" -trunclen 75 -label_suffix _75 -fastqout "${id}_1.trunc.fastq" -padlen 0 -stripleft 0 -stripright 0
  usearch -fastx_truncate "${id}_2.clean.fastq" -trunclen 75 -label_suffix _75 -fastqout "${id}_2.trunc.fastq" -padlen 0 -stripleft 0 -stripright 0

  # Convert to FASTA
  seqtk seq -A "${id}_1.trunc.fastq" > "${id}_forward.fasta"

  # DIAMOND blastx at 75% and 85% identity
  for pid in 85 75; do
    diamond blastx -d plant_associated_db -q "${id}_forward.fasta" \
      -o "${id}_blast_hits_plants_${pid}.txt" --id $pid -e 1e-5 \
      --threads 6 --outfmt 6 --sensitive --query-cover 50

    diamond blastx -d yafL_db -q "${id}_forward.fasta" \
      -o "${id}_blast_hits_yafL_${pid}.txt" --id $pid -e 1e-5 \
      --threads 6 --outfmt 6 --sensitive --query-cover 50
  done

  # Count hits and normalize
  plant_hits_85=$(wc -l < "${id}_blast_hits_plants_85.txt")
  plant_hits_75=$(wc -l < "${id}_blast_hits_plants_75.txt")
  yafL_hits_85=$(wc -l < "${id}_blast_hits_yafL_85.txt")
  yafL_hits_75=$(wc -l < "${id}_blast_hits_yafL_75.txt")
  total_reads=$(grep -c "^>" "${id}_forward.fasta")

  if [[ $total_reads -gt 0 ]]; then
    plant_per_million_85=$(echo "scale=6; $plant_hits_85 * 1000000 / $total_reads" | bc)
    plant_per_million_75=$(echo "scale=6; $plant_hits_75 * 1000000 / $total_reads" | bc)
    yafL_per_million_85=$(echo "scale=6; $yafL_hits_85 * 1000000 / $total_reads" | bc)
    yafL_per_million_75=$(echo "scale=6; $yafL_hits_75 * 1000000 / $total_reads" | bc)
  else
    plant_per_million_85=0
    plant_per_million_75=0
    yafL_per_million_85=0
    yafL_per_million_75=0
  fi

  # Append results
  echo -e "$id\t$plant_hits_85\t$yafL_hits_85\t$plant_hits_75\t$yafL_hits_75\t$total_reads\t$plant_per_million_85\t$yafL_per_million_85\t$plant_per_million_75\t$yafL_per_million_75" >> results_normalized.txt

  # Cleanup
  rm -f "${id}_1.clean.fastq" "${id}_2.clean.fastq" "${id}_1.trunc.fastq" "${id}_2.trunc.fastq"
  pigz -p 6 "${id}_forward.fasta"
done
~~~

Then, to plot the results (Figure 3c), use R:

~~~
> data <- read.delim("/path/to/results_normalized.txt", header = TRUE) #This file should looks like Table S1 in the manuscript.
> library(tidyr)
> library(ggplot2)
> library(dplyr)
> data_long <- data %>%
   pivot_longer(cols = c(YafLHitsPerMillion_85, YafLHitsPerMillion_75, 
                         PlantHitsPerMillion_85, PlantHitsPerMillion_75),
                names_to = "HitType", values_to = "HitsPerMillion")
> ggplot(data_long, aes(x = bulk_vs_rhizos, y = HitsPerMillion, fill = bulk_vs_rhizos)) +
  geom_boxplot(outlier.size = 0.3, lwd = 0.3, alpha = 0.4, size = 0.5,
               linetype = "dashed", width = 0.5, position = position_dodge(0.9)) +
  geom_jitter(color = "black", shape = 21, size = 1, alpha = 0.1,
              position = position_jitter(0.15)) +
  scale_fill_manual(values = c("#904200", "#aad875")) +
  theme_bw() +
  scale_x_discrete() +
  scale_y_log10() +
  ylab("log10(Hits per Million)") +
  facet_wrap(~ HitType, scales = "free_y")
~~~

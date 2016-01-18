# GD_GroSeq_Analysis
GroSeq Documentation for GenomicsData paper.

# Genomics Data Project
## Analysis Documentation

This README file was generated to document all analysis and thought process for this particular project.

## Metagene / Heatmap generation:

Both metagene / heatmaps will be generated using ngsplot version 2.61. I will need to create figures for:

- Pol II
- KAP1
- Hexim
- Larp7
- Cdk9

Bash code and parameters are as follows (You can alter -L parameter for your preference):

### Pol II

    $ ngs.plot.r -G hg19 -R tss -C '/home/carlos/workspace/GenomicsData/data/ChIP-Seq/Pol-II-Chip.mapped.rgid.sorted.filtered.dups_remove.bam' -O pol_ngsplot -F protein_coding -D refseq -L 5000 -P 0 -CO black -GO max -T "Pol II"

### KAP1

    $ ngs.plot.r -G hg19 -R tss -C '/home/carlos/workspace/GenomicsData/data/ChIP-Seq/KAP1/KAP1-Chip.mapped.rgid.sorted.filtered.dups_remove.bam' -O kap1_ngsplot -F protein_coding -D refseq -L 5000 -P 0 -CO purple -GO max -T "KAP1"

### Hexim

    $ ngs.plot.r -G hg19 -R tss -C '/home/carlos/workspace/GenomicsData/data/ChIP-Seq/Hexim/Hexim-Chip.mapped.rgid.sorted.filtered.dups_remove.bam' -O hexim_ngsplot -F protein_coding -D refseq -L 5000 -P 0 -CO red -GO max -T "Hexim"

### Larp7

    $ ngs.plot.r -G hg19 -R tss -C '/home/carlos/workspace/GenomicsData/data/ChIP-Seq/Larp7/Larp7-Chip.mapped.rgid.sorted.filtered.dups_remove.bam' -O larp7_ngsplot -F protein_coding -D refseq -L 5000 -P 0 -CO blue -GO max -T "Larp7"

### Cdk9

    $ ngs.plot.r -G hg19 -R tss -C '/home/carlos/workspace/GenomicsData/data/ChIP-Seq/Cdk9/Cdk9-Chip.mapped.rgid.sorted.filtered.dups_remove.bam' -O cdk9_ngsplot -F protein_coding -D refseq -L 5000 -P 0 -CO darkgreen -GO max -T "Cdk9"

### H3K4me3

    $ ngs.plot.r -G hg19 -R tss -C '/home/carlos/workspace/GenomicsData/data/ChIP-Seq/H3K4me3/H3K4me3.mapped.rgid.sorted.filtered.dups_remove.bam' -O h3k4me3_ngsplot -F protein_coding -D refseq -L 1000 -P 0 -CO orange -GO max -T "H3K4me3"

### H3K27Ac

    $ ngs.plot.r -G hg19 -R tss -C '/home/carlos/workspace/GenomicsData/data/ChIP-Seq/H3K27Ac/H3K27.mapped.rgid.sorted.filtered.dups_remove.bam' -O h3k27ac_ngsplot -F protein_coding -D refseq -L 1000 -P 0 -CO yellow -GO max -T "H3K27Ac"

## Co-occupancy Analysis

Original analysis was done by Jonathan Reeder using a custom made script python script using python 2.7.6. New pairwise co-occupancy analysis requires that I use a version 2.24 of bedtools since MySQL database is not up and thus makes the analysis difficult. Numbers should be similar as previously shown in other trials.

Example of pairwise:

    intersectBed -a Hexim -b Larp7 -u | wc -l

Example of three-way:

   intersectBed -a Hexim -b Larp7 -u | intersectBed -a - -b Cdk9 -u | wc -l

Example of four-way:

   intersectBed -a Hexim -b Larp7 -u | intersectBed -a - -b Cdk9 -u | intersectBed -a - -b KAP1 -u | wc -l

## GRO-SEQ

- Fastq files were taken as documented in the README.md file in the GroSeq folder.

- The three fastq files were then merged together using the cat function:
    
    $ cat '/home/carlos/workspace/GenomicsData/data/GRO-Seq/fastq/SRR828695/SRR828695.fastq' '/home/carlos/workspace/GenomicsData/data/GRO-Seq/fastq/SRR828696/SRR828696.fastq' '/home/carlos/workspace/GenomicsData/data/GRO-Seq/fastq/SRR828729/SRR828729.fastq' > groseq_1_2_3_merged.fastq

- We then mapped the merged fastq file using bowtie v1.1.2 allowing for 2 mismatches:

    $ bowtie ~/libraries/bowtie-1.1.2/indexes/hg19 -n 2 -S -p 8 -q '/home/carlos/workspace/GenomicsData/data/GRO-Seq/fastq/groseq_1_2_3_merged.fastq' > groseq_1_2_3.sam

- The mapped SAM file was then converted to sorted BAM using samtools v1.3:

    $ samtools view -b '/home/carlos/workspace/GenomicsData/analysis/bowtie/groseq_1_2_3.sam' -o groseq_1_2_3.bam

    $ samtools sort -@ 8 '/home/carlos/workspace/GenomicsData/analysis/samtools/groseq_1_2_3.bam' > groseq_1_2_3_sorted.bam

- Sorted BAM file is used as input for bedtools to obtain coverage files (1 for positive and 1 for negative strands):

    $ genomeCoverageBed -bga -5 -strand - -ibam <bamFile> > coverage_neg.bedGraph
    $ genomeCoverageBed -bga -5 -strand + -ibam groseq_1_2_3_sorted.bam > coverage_pos.bedGraph

- With coverage files we can run 'generateMatrix.r' script obtained from Julia of NetSeq paper to get anti-sense / sense read counts. (You simply revert the _neg and _pos nomenclature of the bedgraph files).

### Sense
    $ Rscript '/home/carlos/workspace/GenomicsData/analysis/scripts/generateMatrix.r' coverage_ '/home/carlos/workspace/GenomicsData/data/TSS.txt' groseq_SE.txt

### Anti-sense
    $ Rscript '/home/carlos/workspace/GenomicsData/analysis/scripts/generateMatrix.r' coverage_ '/home/carlos/workspace/GenomicsData/data/TSS.txt' groseq_AS.txt

**NOTES** THE BEDGRAPH FILES _neg _pos NOMENCLATURE IS FOR AS (antisense) READS.

- Use obtained matrix files to run generation of heatmaps and metagene plots as described in Net-Seq paper using script obtained from Julia as above:
    $ Rscript '/home/carlos/workspace/GenomicsData/analysis/scripts/heatmap.r' groseq '/home/carlos/workspace/GenomicsData/analysis/bedgraphtomatrix/' '/home/carlos/workspace/GenomicsData/data/TSS.txt' 

## Divergent Transcripts Analysis
I am attempting to cluster Pol II data to Divergent transcripts in genes. 
## NOTES
Parsing GTF

    awk '{if($3=="genes" && $20=="\"protein_coding\";"){print $0}}' gencode.gtf # 20,345 genes
    awk '$5-$4 >= 2000' file # 18,805 genes
    awk '{gsub(/\"|\;/,"")}1' file

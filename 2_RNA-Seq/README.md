# Basic RNA-Seq processing workflow

## RNA sequencing (RNA-Seq)

&nbsp;&nbsp;&nbsp;&nbsp;RNA sequencing consists on isolating RNA from samples and fragmenting these molecules so that they are fit for [Illumina sequencing](https://www.illumina.com/science/technology/next-generation-sequencing/sequencing-technology.html). Then, these molecules are reversebly transcribed into complementary DNA (cDNA), which is more stable than RNA. cDNA is prepared for sequencing by ligating adapters. Finally, the samples are ready to be sequenced.

![RNA sequencing steps. Source: https://microbenotes.com/rna-sequencing-principle-steps-types-uses/](https://microbenotes.com/wp-content/uploads/2022/07/RNA-Sequencing.jpg)

&nbsp;&nbsp;&nbsp;&nbsp;Sequencing results in FASTQ files, examplified below:

```
@SEQ_ID
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
```

&nbsp;&nbsp;&nbsp;&nbsp;The first is the sequence name (header). The second line is the nucleotide sequence itself. The third line is, optionally, the header again. The last line is consists of sequencing quality values according to the [Phred scale](https://en.wikipedia.org/wiki/Phred_quality_score) for each nucleotide in line 2.

&nbsp;&nbsp;&nbsp;&nbsp;In the case of [paired end sequencing](https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html#:~:text=Paired%2Dend%20DNA%20sequencing%20reads,insertions%2C%20deletions%2C%20and%20inversions.), there are two types of FASTQ files: R1 (forward) and R2 (reverse).

## Workflow parameters

```
params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"

println "reads: $params.reads"
```
# RNA-SEQ Differential Gene Expression Analysis of COVID-19 Infected Pancreas Cells

Next Generation Sequencing technologies promote research in genome wide expression dataand provide high resolution and precise measurements of transcript levels to study gene expression. RNA-Sequencing has become one of the main choices to measure expression levels.
It has many applications such as 'de novo' transcriptome assembly, study of methylation patterns, examining single nucleotide variants, study of alternative splicing etc.

Differential Gene Expression is one of the most common techniques used in downstream analysis of RNA-Seq data. It helps us identify which genes may be under/over-expressed in specific conditions when compared to reference/normal conditions. And these genes whose expression are significantly altered could be the ones contributing to the biological phenomea or disease that we are studying.

In this tutorial, we will be going over a general workflow for Differential Gene Expression using RNA Sequencing data. While these are the most common steps involved in DGE analysis, there are a few things you must keep in mind:
* There are many tools available for each step of DGE analysis. You can choose any of them and this choice is heavily influenced by type of data, organism, sample number, statiscal   requirements as well as goals of research/ analysis.
* RNA Sequencing data analysis is not limited to identifying differential genes alone; it has many applications and these require separate workflows.
* The tools used in this tutorial help analyse the chose dataset but this may not be the same for you. You may have to choose separate tools for your RNA-Seq data.
* When working with large number of samples, it is advisable to work on a high performance cluster or utilise cloud services.

Let's begin the tutorial with a basic introduction to RNA Sequencing.

## Table of Contents

- [Introduction](#intro)
     - [Introduction to RNA sequencing](#general_intro)
     - [Background of the dataset](#dataset_intro)
     - [Installation](#installation)
     - [Setting up Working Directory](#work_directory)
     - [Downloading Data using SRA toolkit](#download_data)
- [Differential Gene Expression Analysis from RNA-Seq Data](#workflow)
     - [Quality Control of RNA Seq Reads](#quality_control)
     - [Preprocessing of RNA Seq Reads](#preprocessing)
     - [Alignment of RNA Seq Reads to Reference Genome](#alignment)
     - [Transcriptome assembly of Aligned Reads](#assembly)
     - [Reestimating Gene Counts Matrix](#gene_counts)
     - [Differential Gene Expression Analysis](#diff_gene)
- [Citations](#citations_list)
- [License](#license_name)

## Introduction <a name="intro"></a>

### What is RNA Sequencing? <a name="general_intro"></a>

RNA Sequencing is a bundle of experimental (wet-lab) and computational techniques that help determine the identity and abundance of RNA sequences in biological samples. 
It provides access to the transcriptome – complete set of transcripts present in a cell at a given developmental state or physiological condition.

It has several advantages over earlier technologies such as microarrays: high-throughput, increased sensitivity that provides minute details of transcriptional features and detection of novel transcripts, gene models and small noncoding RNA species.

In the wet-lab procedures, RNA is isolated from freshly dissected or frozen cells or tissue samples. Quality control practices are integrated to check for degradation, purity and quantity.

The quality-checked total RNA is used for library preparation where RNAs in the sample are converted to cDNA (complementary DNA) by reverse transcription. cDNAs are more stable and amenable to the sequencing technology. The libraries are then to sequencing facilities where different sequencing protocols are used depending on the platform.

The result of sequencing are raw reads which are further analysed with bioinformatic and computational techniques. In the next section, we will go over the different steps involved in the computational aspect of RNA-Seq analysis. 

### What are the general steps involved in RNA Sequencing Data Analysis?

Millions of reads are obtained from the RNA sequencing experiments which are analysed with the help of computational approaches and statistics. These techniques assist in gene and splice variant discovery, differential expression analysis and detection of fusion genes, variants.
The following image gives you a brief idea of the computational approaches involved in RNA-Seq data analysis.


<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/rna_seq_workflow.png" width="400" height=600 alt="RNA Sequencing Workflow image"/>
</p>

<p align="center">
     <b>Differential Gene Expression in RNA Sequencing - Workflow </b>
</p>

On the left are the basic steps involved in RNA Sequencing Data Analysis and on the right are the tools/software that are used in these procedures. 
Not all analysis, ends with Differential Expression Analysis. Sometimes pre-processed reads can be directly used for de-novo transcriptome assembly followed by annotation.

### What tools are will be using in this RNA Sequencing Analysis - Tutorial?

For this tutorial, we will be using the following tools
- Quality Control : FASTQC, RNA-SeQC
- Pre-processing: Trimmomatic
- Read Alignment: HISAT2
- Transcriptome Assembly: StringTie
- Differential Analysis: DESeq2

Please download these tools and install them in your system. Before proceeding with the tutorial, please install Ubuntu if you are using Windows OS. We will be using Bash/Shell and RStudio in this tutorial.

### Which dataset will we use for DE gene analysis from RNA Seq data? <a name="dataset_intro"></a>

The dataset we will be using in this tutorial is GSE159717 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159717) which can be downloaded from Gene Expression Omnibus (GEO) (https://www.ncbi.nlm.nih.gov/geo/)
This dataset consists of samples involving SARS-CoV-2 infection in human pancreatic islet cells. The datasets consists of 6 RNA-seq samples:
-	2 control samples of human pancreatic islets
-	2 samples of human pancreatic islets with SARS-CoV-2 infection
-	2 samples of human pancreatic islets with SARS-CoV-2 infection but treated with Remdesvir
              
We will need the FASTQ files for the analysis. This can be downloaded in 2 ways:

- From NCBI’s Sequence Read Archive Using the SRA toolkit: 
You can download the SRA toolkit (https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) and install it. 
SRA toolkit can be run via command line. You can download each FASTQ file using SRR ID and the commands prefetch and fasterq-dump. We will cover this in the tutorial.

- From ENA browser: You can go to the European Nucleotide Archive (https://www.ebi.ac.uk/ena/browser/) and enter the accession number PRJNA670242 on the top right.
Click on the view button and that will lead you to a page giving details about the BioProject PRJNA670242. The table on the page enlists all the 6 FASTQ files that can be downloaded by clicking on the ‘Download All’ button.

### What is the scientific background of the dataset?

The dataset GSE159717 involves RNA Sequencing on samples of human pancreatic islet cells obtained from patients diagnosed with COVID-19.
This study was conducted by Müller and colleagues involving Ulm University Medical Center, Ulm, Germany and multiple institutions. Please check this link (https://pubmed.ncbi.nlm.nih.gov/33536639/) for the full affiliation list.

***Citation**: Müller JA, Groß R, Conzelmann C, Krüger J, Merle U, Steinhart J, Weil T, Koepke L, Bozzo CP, Read C, Fois G, Eiseler T, Gehrmann J, van Vuuren J, Wessbecher IM, Frick M, Costa IG, Breunig M, Grüner B, Peters L, Schuster M, Liebau S, Seufferlein T, Stenger S, Stenzinger A, MacDonald PE, Kirchhoff F, Sparrer KMJ, Walther P, Lickert H, Barth TFE, Wagner M, Münch J, Heller S, Kleger A. **SARS-CoV-2 infects and replicates in cells of the human endocrine and exocrine pancreas.** Nat Metab. 2021 Feb;3(2):149-165. doi: 10.1038/s42255-021-00347-1. Epub 2021 Feb 3. PMID: 33536639.*

COVID was initially considered an exclusive lung disease. However, as the infection spread and cases increased, evidences from clinical and experimental studies show that the virus damages other organs such as kidneys, heart, brain, gastrointestinal and endocrine organs.

In this study, the researchers focused on SARS-Cov-2 infection in the pancreas. The pancreas is an organ located in the abdomen that aids in conversion of food to energy. It has:
 * Exocrine function: It releases digestive enzymes in the small intestine to help break down the food into fats, carbohydrates and proteins.
 * Endocrine function: It releases hormones such as insulin and glucagon into the bloodstream that helps maintain optimal blood sugar levels.

Multiple studies showed a connection between COVID infection and poor pancreatic function. 
 * Pre-existing diabetes increased the risk of developing SARS-Cov-2 infection; it required intensive treatment and was associated with increasing mortality.
 * SARS-Cov-2 infection affected the exocrine function of the pancreas leading to pancreatitis (severe inflammation of pancreas), pancreatic enlargement and abnormal levels of    digestive enzymes in the patients.
 * Increased blood sugar levels was observed in patients with type2 diabetes and SARS-Cov-2 infection
 * Ketoacidosis (increase in ketone levels due to low insulin levels which makes the blood acidic) was observed in diabetic and non-diabetic patients with SARS-Cov-2 infection
 * New onset of type I diabetes in the absence of autoantibodies following SARS-Cov-2 infection was observed. 

SARS-Cov-2 might trigger beta-cell injury (cells responsible for production of insulin) either by disrupting immune function or by directly interrupting beta cell function.
The researchers pointed out that several cellular factors or proteins expressed on target cells could be responsible for entry of SARS-Cov-2 in pancreatic cells that lead to subsequent destruction. These factors include:
* **angiotensin-converting enzyme 2 (ACE2)**: protein that affects blood sugar levels and blood pressure and it also serves as a receptor for SARS-Cov-2 viral entry.
* **transmembrane serine protease 2 (TMPRSS2)** : protein that is expressed on cells of multiple organs that assists in breakdown of proteins which has multiple implications in body functions such as immune function and allergies, development of prostate cancer and influencing exocrine function of the pancreas. TMPRSS2 also assists in viral entry and spread.

To examine the possible mechanism of viral entry and eventual disturbance of pancreatic function, the researchers conducted the multiple experimental studies; however, we will be focussing on RNA sequencing analysis on uninfected and infected (with or without remdesivir), cultured human islets from two donors utilised in the above mentioned study. Expression profiling libraries were sequenced on a HiSeq 4000 instrument (Illumina) in 50-bp, single-end mode.

---------------------------------------------------------------------------------
### Installation of Packages <a name="installation"></a>

This tutorial assumes that you have the basic knowledge of Linux/Shell scripting and R programming. To implement this tutorial, please ensure that you have the following installed in your system:

- If you are using a Windows system, install Ubuntu
- R 4.0.0
- Python 3
- RStudio
- FastQC
- HISAT2
- StringTie
- SAMtools

The R packages that will be required later on in the analysis can be installed using RStudio and the installation procedure has been described further on in the tutorial.

--------------------------------------------------------------------------------------

### Setup your working directory <a name="work_directory"></a>

In the image above, it is pretty obvious that RNA-SEQ DEG analysis consists of multiple steps involving different input and output inputs. Therefore it is very important to organize your data into separate folders, so that you don't get lost while working. :)

Here is how you can organize your working directory structure:

``` bash
── rna_seq_dge_analysis/
  │   └── alignment/                    <- Data generated during alignment steps
  │       ├── 1_index/                  <- Folder to store the indexed genome files from HISAT2
  │       ├── 2_output/                 <- Alignment files generated from HISAT2 (.SAM)
  │       ├── 3_input_reads/            <- Trimmed and QC checked input reads for each sample for genome alignment (.FASTQ)
  │   
  │    └── assembly/                     <- Data generated during transcriptome assembly steps
  │       ├── 1_annotation/             <- Folder to store Genome annotation file (.GTF/.GFF)
  │       ├── 2_merged_transcripts/     <- Data generated during stringtie merge_transcript step
  │       ├── 3_output/                 <- Folder to store transcriptome assembly output files (.GTF)
  │       ├── 4_sorted_bam_reads/       <- Folder to store sorted BAM input files for transcriptome assembly (.BAM)
  │
  │    └── data/                         <- Location of input  RNAseq data
  │  
  │    └── qc/                           <- Data generated during quality control and pre-processing steps
  │       ├── 1_fqc_results/            <- Results of FASTQC for each sample
  │       ├── 2_trim_fastqc/            <- Results of FASTQC for every trimmed sample
  │       ├── 3_trim_output/            <- Output files of Trimmed reads for every sample
  │
  │    └── ballgown/                     <- Data generated during Stringtie's quantification step; stores transcript quanitifcation data for every sample
  │  
  │    └── dge_analysis/                 <- Folder to store files generated during DESeq2 analysis      
  │       
   
``` 
--------------------------------------------------------------------------------------------------

### Downloading Input Read Files with SRA Toolkit <a name="download_data"></a>

To begin our RNA-seq data analysis, we need the raw FASTQ reads for each sample. As mentioned previously, you can download these from European Nucleotide Archive with the accession number PRJNA670242.
Else you can use SRA Toolkit. Assuming that you have downloaded and installed SRA Toolkit (https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software), open your bash terminal and add the location/address of sra-toolkit (the place where it is installed) to your path.

the command to add the address to the path:
export PATH=$PATH:$PWD/sratoolkit.2.10.8-ubuntu64/bin

export PATH=$PATH:/place/with/the/file   [General usage]

To avoid doing this everytime, you can add the path Global shell specific configuration files or Per-user shell specific configuration files. You can read more about it here https://linuxize.com/post/how-to-add-directory-to-path-in-linux/

Once that's done, navigate to your working directory.

```bash

cd rna_rna_seq_dge_analysis
cd data
prefetch SRR12852623

```

You should see the following output:

```bash

2022-01-05T05:25:08 prefetch.2.10.8: 1) Downloading 'SRR12852623'...
2022-01-05T05:25:08 prefetch.2.10.8:  Downloading via HTTPS...
2022-01-05T05:35:51 prefetch.2.10.8:  HTTPS download succeed
2022-01-05T05:35:53 prefetch.2.10.8:  'SRR12852623' is valid
2022-01-05T05:35:53 prefetch.2.10.8: 1) 'SRR12852623' was downloaded successfully
2022-01-05T05:35:53 prefetch.2.10.8: 'SRR12852623' has 0 unresolved dependencies

```

Next run the following command:

```bash

fastq-dump SRR12852623

```

You should see the following output:

```bash

Read 27218187 spots for SRR12852623
Written 27218187 spots for SRR12852623

```

You should run the same commands for each of the samples. I have enlisted the commands below, but you have to wait after running each command for the output and then enter the next.

```bash

prefetch SRR12852624
fastq-dump SRR12852624

prefetch SRR12852625
fastq-dump SRR12852625

prefetch SRR12852626
fastq-dump SRR12852626

prefetch SRR12852627
fastq-dump SRR12852627

prefetch SRR12852628
fastq-dump SRR12852628

```
---------------------------------------------------------------------------------

## Differential Gene Expression Analysis from RNA Sequencing Data - Workflow <a name="workflow"></a>

Now, we will begin our data analysis. One important thing to note is that this is an exploratory data analysis; the number of samples/ replicates in this dataset are small. They are great to find some new leads on how COVID-19 infection affects pancreatic cell functions and diabetes onset but we definitely need larger number of samples to make definitive conclusions.
Nevertheless, the aim of this tutorial is to give you a general idea of DGE RNA Seq analysis and this dataset helps achieve just that.
Lets's begin our analysis with quality control of our raw RNA-Seq reads.

---------------------------------------------------

### 1. Quality Control of RNA-Seq Samples <a name="quality_control"></a>

Quality-related issues generally arise during sequencing or library preparation. These include: low confidence bases, PCR artifacts, sequence-specific bias, sequence contamination, untrimmed adapters, 3’/5’ positional bias.
Running quality checks can help avoid problems that would occur during genome alignment. Few quality-related problems can be corrected. Some cannot be corrected, but being aware of them can help us interpret results cautiously.

FASTQC is a Java program that performs multiple quality checks on tens of millions of reads in a few minutes. It reports and helps visualize information on base content and quality, k-mer content, presence of ambiguous bases, overrepresented sequences, duplicates etc.

Add FASTQC to your path, if not yet added to the configuration file. Navigate to your working directory: rna_seq_dge_analysis and then run this command.

The following command produces a quality report for a single sample:

```bash

fastqc ./data/SRR12852623.fastq -o ./qc/fqc_results

```

You can run the same command (but replacing the name of the FASTQ file) for each sample or you can run this command for all samples at once

```bash

fastqc ./data/*.fastq -o ./qc/fqc_results

```

Let’s look at some quality check reports produced by FASTQC for sample SRR12852623.fastq
<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/fastqc_image1.png" width="800" height=400 alt="FastQC basic statistics image"/>
</p>

<p align="center">
     <b>FastQC Report: Basic Statistics for sample SRR12852623 </b>
</p>

At the left of the image, FASTQC has given judgements (pass, warn, fail) on several quality metrics. These judgements are based on general thresholds and poor judgements may not always suggest that the sample has failed quality checks.
For example, the sample fails ‘Sequence Duplication Levels’ but this generally acceptable in RNA-Seq samples.

In the Basic Statistics section, we can see that the sample SRR12852623.fastq has 27218187 sequences with each read length of 50 bases. Also the base quality is given in Sanger/Illumina 1.9 encoding or Phred 33 score.
<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/fastqc_image2.png" width="800" height=400 alt="FastQC per base sequence quality image"/>
</p>

<p align="center">
     <b>FastQC Report: Per base sequence quality for sample SRR12852623 </b>
</p>

This plot shows the range of quality values across all bases at each position in the FASTQ file. Base quality indicates the confidence in the base call or how correctly the sequencer has identified the base at a given position in the given sequence.

These scores are expressed in Phred scale which is given as Q = -10 log10 P where P is the probability that the base is wrong. The scores range from 0 to 40. In the FASTQ files, they are encoded as ASCII characters instead of numbers to save space.

Phred 33 score refers to the encoding where the 33rd ASCII character is used as 0. Phred 64 score is used in old Illumina software the 64th ASCII character is used as 0.

In the Base quality report, at each position a BoxWhisker type plot is drawn. Here the upper and lower whiskers indicate 10% and 90% points while blue line indicates mean quality and the red line indicates median quality.

The background of the graph divides the y axis into very good quality calls (green), calls of reasonable quality (orange), and calls of poor quality (red). 
For the sample SRR12852623, we can see uniformly high-quality scores across all positions in the sequences.
<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/fastqc_image3.png" width="800" height=400 alt="FastQC per sequence quality image"/>
</p>

<p align="center">
     <b>FastQC Report: Per sequence quality score for sample SRR12852623 </b>
</p>

The per sequence quality score plot examines average quality score over the full length of the read for a subset of sequences. Here majority of the reads should have a high average quality score with no large bumps at the lower quality values.

In this plot, we can see that average quality score per read is 40 for a subset of 1.4 x 10<sup>7</sup> reads.
<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/fastqc_image4.png" width="800" height=400 alt="FastQC per base sequence content image"/>
</p>

<p align="center">
     <b>FastQC Report: Per base sequence content plot for sample SRR12852623 </b>
</p>


The per base sequence plot reports the percent of bases called at each position across all reads in the file. The ‘Per Base Sequence Plot’ shows a Fail in the report for sample SRR12852623.fastq 

This is acceptable for RNA-Sequencing; a non-uniform distribution for the first 10-15 bases of the read. This is due to random hexamer priming that is utilised in RNA-Sequencing library preparation.

Random hexamer primers are a mix of oligonucleotides representing hexamer sequences that are attached to the single stranded RNA for extension by reverse transcription. 

While this process is intended to be random, multiple studies have shown that random hexamer priming introduces a bias in the nucleotide composition at the start of the reads.
This affects the expression estimates of genes and isoforms and the also the resulting sequence coverage of the transcripts is not uniform.

In most cases, this does not affect downstream analyses.
<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/fastqc_image5.png" width="800" height=400 alt="FastQC per sequence GC content image"/>
</p>

<p align="center">
     <b>FastQC Report: Per sequence GC content plot for sample SRR12852623 </b>
</p>

The ‘Per Sequence GC content’ average GC content over all sequences (indicated as red line) and compares it modelled normal distribution of GC content. 

In a normal random library, a roughly normal GC content would be observed that corresponds to overall GC content of the genome of that particular organism. Since the original GC content is not known, a model of the GC content is developed based on the reference data (indicated as a blue line).

In this plot, we observe that both the distributions are almost similar. In case, the red distribution would be unusually different from the blue one, this could mean that the library is contaminated with a genome of another organism (could be observed by broad peaks) or there are other kinds of bias (in case of over-represented sequences, you would see sharp peaks). 
<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/fastqc_image6.png" width="800" height=400 alt="FastQC sequence duplication levels image"/>
</p>

<p align="center">
     <b>FastQC Report: Per sequence duplication levels for sample SRR12852623 </b>
</p>

The ‘Sequence Duplication Levels’ plot shows the relative number of sequences with different degrees of duplication.

This module analyses only first 100,000 sequences in each file. Each sequence is tracked to the end of the file to give a representative count of the overall duplication level.

Reads above 75bp are truncated to help in detecting exact match of the duplicated sequences. Also long sequences tend to have sequencing errors that would make them diverse and underrepresent the duplication levels.

In a diverse library most sequences will occur only once in the final set. A low level of duplication may indicate a very high level of coverage of the target sequence, but a high level of duplication could indicate a bias such as PCR over-amplification or low complexity library due to small amounts of starting materials.

In RNA-Sequencing, duplicates are often a natural consequence of sequencing highly expressed transcripts. Thus it is highly likely that this plot would ‘FAIL’ for RNA-Seq samples.
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/fastqc_image7.png" width="800" height=400 alt="FastQC overrepresented sequences image"/>
</p>

<p align="center">
     <b>FastQC Report: Overrepresented sequences for sample SRR12852623 </b>
</p>

The ‘Overrepresented Sequence Plot’ lists all of the sequence which make up more than 0.1% of the total. 

Theoretically, a normal-high throughput library would have a diverse set of sequences; no individual sequence would account for a high fraction of the whole. However, if the sample does contain overrepresented sequences, it could mean that the plot is highly biologically significant or the library is contaminated or it has a bias.

The overrepresented sequences could be vector or adapter sequences. One can BLAST the sequence to determine the identity.

In RNA-Seq, some transcripts may be so abundant that can be listed as overrepresented sequences. Also in case of small libraries, where sequences are not subjected to random fragmentation, one sequence may account for a huge proportion of total reads.

----------------------------------------------------------------

### 2. Preprocessing of RNA-Seq samples <a name="preprocessing"></a>

After conducting quality checks, multiple pre-processing steps can be conducted to mitigate some of the quality problems that arose during the experimental setup. 
This assists in better alignment of the reads to the genome.

These steps include:
* **Filtering**: You can filter reads based on their quality. Average read quality is calculated and if it falls below the user-defined threshold that read is dropped from further analysis. In case of Paired End reads, if the mean quality of either of the reads drops below the given threshold, the read pair is dropped from further analysis.
* **Trimming**: Instead of dropping entire reads or read pairs from further analysis, you can trim low quality bases from a given read. There are several ways to trim a read:
   * *Trimming from either end*: Starting from either 3’ or 5’ end, bases are trimmed if their quality falls below the user-defined threshold. If the trimmed reads fall below a user-defined length, they can be filtered.
   * *Trimming using a sliding window approach*: Instead of examining the quality of one base at a time, you can define a length of a search window and examine the mean quality of bases in that window. If the mean quality is above the given threshold, the window slides further to examine the next set of bases.</br>
Sliding the window from the 5′ end keeps the beginning of the read until the quality falls below the defined threshold, while sliding from the 3′ end cuts until it reaches a window with good enough quality. 
The window size is an essential parameter that needs to be tuned; very small window size may lead to a stringent check and lead to loss of reads.
   * *Trimming by sum method*: This is also known as BWA quality trimming. The read is scanned from the 3’ end, the quality of each base is compared to the user defined threshold and the difference is summed up. The read is trimmed at the point where the ‘summed up difference’ is the highest.
* **Removal of Ambiguous Bases**: If a base is not identified during sequencing, it is indicated as a N. Higher number of Ns in the read should be removed to avoid false mapping and incorrect transcriptome assemblies.
* **Removal of Adapters**: Adapters are short, known sequences of oligonucleotides that are used to extract or fish out DNA sequences of interest. Other tags such as primers and multiplexing identifiers may also be attached to the reads.
These need to be removed prior to further analysis. However, removal of adapter content can present several hurdles.</br>
Like any other part of the read, these tags can undergo sequencing errors like mismatches, indels and ambiguous bases. When sequencing small-RNAs, reads can run into a ‘read-through’ situation where the reads extend into the adapter and 3’ end adapter can be partial.</br>
Trimming tools overcome these challenges by examining the reads for known adapter content, aligning the portions of the reads that partially match these adapter sequences and then trimming the matched portions.
* **Examining Read Length**: Read length distribution gives us an idea of how useful the reads are for further steps such as genome alignment, transcriptome assembly and detection of splice isoforms. Very short reads, resulting from trimming and adapter removal, map unambiguously to the genome and hence can be removed from further analysis.</br>

The above-mentioned steps are the most common procedures utilised for pre-processing data prior to alignment. Depending on your data quality and the QC plots, you may need to take additional steps such as cautious removal of duplicates, examining reads for possible contaminants and removing related sequences, removal of low-complexity sequences and poly A/T tails etc.

Several tools are available for pre-processing the data: Trimmomatic, PRINSEQ, Cutadapt, TrimGalore, FastX, TagCleaner 

In this tutorial, we will be using Trimmomatic. We have some really good quality reads in our dataset, so it is quite likely that trimming may not have a significant effect on the read quality.
However we will trim and remove adapters from our samples and then examine the QC reports once again.

Assuming you have downloaded and installed Trimmomatic. (http://www.usadellab.org/cms/?page=trimmomatic) Add Trimmomatic to your PATH or configuration file.

```bash

TrimmomaticSE -phred33 ./data/SRR12852623.fastq  ./qc/trim_output/SRR12852623_trimmed.fastq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:40

```
In this command, the options are:
* phred 33 : indicates that base quality is given Phred33 format
* ILLUMINACLIP: TruSeq3-SE.fa:2:30:10 : it points to the FASTA file containing adapter sequences that are used in Illumina MiSeq and HiSeq protocols. It allows two mismatches in the seed, palindrome clip threshold is 30, simple clip threshold is 10
* SLIDINGWINDOW:5:20 : indicates that a sliding window of 5 bases is to be used and reads are to be trimmed if quality of the window drops below 20
* MINLEN:40 : indicates that minimum read length should be 40 after trimming

You have to run the above mentioned commmand for each of your samples: SRR12852624, SRR12852625 ...etc; just replace the name of the sample. Or you could run the following bash script to process all samples together

```bash
#!/bin/bash

SAMPLES="SRR12852623 SRR12852624 SRR12852625 SRR12852626 SRR12852627 SRR12852628"

for SAMPLE in $SAMPLES; do
    TrimmomaticSE -phred33 ./data/${SAMPLE}.fastq ./qc/trim_results/${SAMPLE}_trimmed.fastq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:40

done

```
This script is available as samples_trim.sh in <u> <b> the scripts folder in the repository. </b> </u>

After trimming our samples, we can run quality control checks to revaluate them; however since we started with really high quality reads, we most likely won't see any changes. 

```bash

fastqc ./qc/trim_output/*.fastq -o ./qc/trim_fastqc

```

Open the fastqc html for each sample. After running FASTQC again on the trimmed reads, we observe no significant changes in the reports except for ‘Sequence Length distribution’ plot. After trimming, reads are likely to be of different lengths.

-----------------------------------------------------------------

### 3. Genome Alignment of Preprocessed RNA-Seq reads <a name="alignment"></a>

Genome alignment involves mapping reads to a reference genome. This helps us identify where the reads originate from and how similar they are to the reference genome.
This data can be used to discover new genes and transcripts and for quantifying expression. In the absence of a reference genome, reads can be mapped to a transcriptome (collection of all RNA transcripts present in the cell of a given organism at a given time).

However, there are a couple of challenges when mapping reads to a genome:
* There are millions of reads and they are short in length. Genomes are large and do contain sequences like repeats or pseudogenes which make it difficult to map each read to a unique position.
* Additionally, there are mismatches and indels caused by genomic variation and sequencing errors that need to be dealt with.
* Many organisms have introns (non-coding areas of an RNA transcript) in their genes, so the reads may align non-continuously to the reference genome.

Multiple alignment programs have been developed to overcome these challenges. Here we will be using HISAT2 (Hierarchical Indexing For Spliced Alignment of Transcripts).

It is a spliced aligner – used for aligning spliced reads. Introns are non-coding sequences of mRNA transcripts while exons are coding sequences. A single gene can have multiple exons that again can be arranged in multiple different combinations resulting in diverse mRNA transcripts and different proteins. This is called alternative splicing.

HISAT is based on Burrows Wheeler Transform (BWT) which is a data compression algorithm. It transforms data in a way that it is amenable to compression and is useful for data containing lots of repeats (sequence information, in our case).

HISAT2 extends the idea of BWT for graphs. It creates a graph-based index of the reference genome. Graph BWT involves finding matching paths in the graph that are equivalent to certain sections of the reads.

HISAT2 builds a whole genome global index and multiple small local indexes to make spliced alignment possible. With the human genome, HISAT2 builds one global index and 48000 local indexes
Additionally, the graph-based index captures a wide representation of genetic variants (SNPs) with very low memory requirements. 

The basic approach of the aligners is ‘seed and extend’ which involves:
1. identifying segments of reads of defined lengths (seeds) that precisely map to a given location in the genome. Seed matches can be exact or tolerate mismatches.
2. extend the reads in both directions to map rest of the read or maximum mappable length

HISAT2 maps longer part of the read that maps to the genome contiguously using the global index. Once this is mapped, it helps to identify the relevant local index.  HISAT2 can usually align the remaining part of the read within a single local index rather than searching across the whole genome.

You can download and install HISAT2 from http://daehwankimlab.github.io/hisat2/download/ . Additionally, download the reference genome index (H.sapiens GRCh38 genome) given by HISAT2 developers so as to save time.

Unzip and save the downloaded ‘genome’ index in the index folder under the alignment directory. 

Navigate to your working directory rna_seq_dge_analysis. Copy the trimmed reads to the input_reads folder in the alignment directory using the following command: 

```bash

cp ./qc/trim_output/*.fastq ./alignment/input_reads

```

Assuming you have downloaded and installed HISAT2. Add HISAT2 to your path, if not added to the configuration file. Now run HISAT2 using the following command:

```bash

hisat2 –p 4 -f –x ./alignment/index/genome -q -U ./alignment/input_reads/SRR12852623_trimmed.fastq –S ./alignment/output/SRR12852623_aligned.sam

```

You have to run this command for each sample by changing the sample names. Alternatively, you can use the samples_alignment.sh script present in **the scripts folder in the repository.** In this command -p refers to the number of processors; here I have used 4 since I have access to it, you can change this number depending on your PC configuration.

---------------------------------------------------------------------

### 4. Transcriptome Assembly with Aligned RNA-Seq reads <a name="assembly"></a>

Transcriptome assembly is performed to obtain full-length transcripts based on the sequence reads. DNA sequence of one or more genes is transcribed into RNA and this is referred to as RNA transcript. A mature RNA transcript comprises of a combination of exons (protein coding sequences).

Transcriptome assembly is different from genome assembly. Genome assembly generally has uniform read coverage , making exceptions for sequencing and library preparation biases. A deviation from uniform coverage could indicate presence of repeats.

But, in transcriptome assembly or RNA-seq data, the number of transcripts between genes or various isoforms of same gene could vary by several magnitudes and are expressed at different levels. This difference in expression levels could contribute to the non-uniformity observed in read coverage data.

There are two ways of performing transcriptome assembly:
1. Reference guided assembly:  Reads are first mapped on the genome and the assembly task consists of solving which mapped reads correspond to which transcripts.
2. De-Novo assembly: In the absence of a reference genome, the assembly is based on utilizing sequence similarity between the RNA-seq reads. 

We will be utilising the first step- Reference guided assembly. We will use StringTie to perform transcriptome assembly.

StringTie2 is a fast assembler of RNA-Seq alignments into potential transcripts. As a reference guided assembler, it takes advantage of an existing genome to which the reads are aligned. It then builds splice graphs based on these alignments and uses the graphs to construct individual transcripts.

It also offers optional de novo assembly step to assemble and quantitate full-length transcripts representing multiple splice variants for each gene locus.

You can download StringTie from https://ccb.jhu.edu/software/stringtie/index.shtml. Add it to your PATH.

Before we begin transcriptome assembly with StringTie2, we need to preprocess our SAM files obtained in the previous step. 

#### 4a. Converting SAM files to BAM files

We have to convert SAM files to BAM files and sort them by coordinate before using them further. We will use SAMtools for this step. 

SAM format is Sequence Alignment Map format that stores biological sequences that have been aligned to a reference sequence. The binary equivalent of SAM is BAM (Binary Alignment Map) which stores the same data in a binary representation to save space.

SAMtools is a package that helps manipulate SAM format: such as SAM/BAM conversion, sorting, indexing, or merging. Please download and install SAMtools from http://www.htslib.org/  and add it to your path.

We will convert the SAM files to BAM and sort the alignments by chromosomal coordinates using the following command:

```bash

samtools sort -@ 4 -o ./assembly/sorted_bam_reads/SRR12852623_sorted.bam ./alignment/output/SRR12852623_aligned.sam

```

Run this for each sample by replacing the sample name or you can use the script samples_sortbam.sh from **the scripts folder in the repository.** Here -@ refers to the number of processers and you can change that as per your PC configuration.


#### 4b. Reference-Guided Transcript Assembly

To begin reference-guided transcript assembly, we first need the reference human genome GTF file. GTF refers to Gene Transfer Format. It contains information about gene structure such as chromosome ID, annotations from public database, features like coding sequence, intron, exon, start codon, stop codon, 5’ UTR, 3’ UTR etc.

You can download the GTF file for human genome from Ensembl. (https://asia.ensembl.org/info/data/ftp/index.html) Download the file Homo_sapiens.GRCh38.104.chr.gtf.gz , extract and save it in the annotation folder in assembly directory. Rename it as ‘Homo_sapiens_chr.gtf’ for the sake of convenience.

Navigate to your working directory rna_seq_dge_analysis. Now, we can assemble the transcripts for each sample using the following command:

```bash

stringtie -p 4 -G ./assembly/annotation/Homo_sapiens_chr.gtf -o ./assembly/output/SRR12852623.gtf -l SRR12852623 ./assembly/sorted_bam_reads/SRR12852623_sorted.bam

```

You can run the same command for each sample or use the script samples_transcriptassembly1.sh from **the scripts folder in the repository.**


#### 4c. Merging Transcripts

Transcripts assembled for all samples can be merged together for further analysis using the following command. This is extremely beneficial when you have hundreds or thousands of samples.

Navigate to your working directory rna-seq-dge-analysis and run the following command:

```bash

stringtie --merge -p 4 -G ./assembly/annotation/Homo_sapiens_chr.gtf -o ./assembly/merged_transcripts/stringtie_merged.gtf ./assembly/mergelist.txt

```

Again, you can rest the parameter 'p' to the number of processors in your PC. Also mergelist.txt file is a notepad file containing paths to every GTF file for each sample. It is provided in **the tutorial_data folder in the repository.**

The output of this command id stringtie_merged.gtf file.


#### 4d. Re-evaluating StringTie transcript abundance estimates

Finally, we redo Stringtie abundance calculations using the newly merged transcripts. This needs to be done to generate count files for each sample.

Assuming you are in your working directory rna_seq_dge_analysis, run this command for a single sample:

```bash

stringtie -e -p 4 -B -G ./assembly/merged_transcript/stringtie_merged.gtf -o ./quantification/SRR12852623/SRR12852623_requant.gtf ./assembly/sorted_bam_reads/SRR12852623_sorted.bam

```

You can run this for each your samples or you can use the script samples_transcriptassembly2.sh present in **the scripts folder of the repository.**

#### 4e. Generating Count Tables for Genes and Transcript estimates

Now we have the quantification folder ready for differential expression analysis. Since we will use DESeq2 for DGE analysis, we will use the prepDE.py3 (python 3 script) to generate count tables for genes and transcripts.

This script is provided by the Stringtie authors. You can download the prepDE.py3 file from https://ccb.jhu.edu/software/stringtie/dl/prepDE.py3

Next create a text file that contains the SAMPLE ID followed by <path to SAMPLE ID> and name it samples_lst.txt You can find this file in **the tutorial_data folder in the repository.**

Install python3 for Linux (if its not already installed). You can follow this guide https://docs.python-guide.org/starting/install3/linux/

Save both prepDE.py3 and samples_lst.txt file in your working directory rna_seq_dge_analysis. Run the following command:

```bash

python3 prepDE.py3 -i samples_lst.txt 

```

This will generate two files gene_count_matrix.csv and transcript_count_matrix.csv  Move these files to dge_analysis folder. 
 
--------------------------------------------------------------------

### 5. Reestimating Gene Counts from Assembled RNA-Seq reads <a name="gene_counts"></a>

We have assembled our RNA-Seq alignment into potential transcripts and we have also generated read count tables at both gene and transcript levels.
The files are gene_counts.csv and transcript_counts.csv.

We ***WOULD*** (actually, we aren't) be using the gene_count_matrix.csv file for differential gene expression analysis. A closer look at the contents of gene_count_matrix.csv file reveals that we have 79,000+ entries.

Most are labelled with ENSEMBL identifiers 'ENSG' but many of them have a unknown identifier or label starting with 'MSTRG'. MSTRG IDs are default names assigned by Stringtie while merging transcript gtfs. 

If we were to use this data for further analysis and eliminate entries that do not have ENSEMBL identifiers we might miss out on valuable gene entries.

Here is where IsoformSwitchAnalyzeR comes to our aid. 

#### What is IsoformSwitchAnalyzeR?

IsoformSwitchAnalyzeR is an R package that is designed to identify isoform switches with the help of statistics from RNA sequencing derived quantification of novel or annotated (known) isoforms.
It accepts input data from tools such as Cufflinks, StringTie, Kallisto and Salmon.

Each gene produces different transcripts or isoforms with the help of Alternative Splicing, alternative transcription start sites (aTSS), and alternative transcription 
termination sites (aTTS). 

Many of these isoforms tend to have different functions. Isoform switching is 'the differential usage of isoforms under different conditions' and since isoforms tend to differ in function, their altered usage can have a significant biological impact.

Different isoform usage is commonly observed in normal biological processes such as cell development, pluripotency and apoptosis but their dysregulated use can be responsible for conditions like cancer. 

IsoformSwitchAnalyzeR helps conduct genome-wide analysis of specific types of alternative splicing and also predicts functional consequences of isoform switches.

You can learn more about it here https://bioconductor.org/packages/release/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html 

#### Using IsoformSwitchAnalyzeR to rescue StringTie Annotation
However, we will not be conducting an Isoform Switch Analysis. IsoformSwitchAnalyzeR has this novel algorithm that helps 'rescue StringTie annotation and extract gene count matrix.' 

The authors of IsoformSwitchAnalyzeR highlight that StringTie provides its own gene ids 'MSTRG.XXXX' for many genes. These are sets of overlapping transcripts and are used instead of known 'reference gene ids' because:
* For all genes, when a 'novel' transcript is identified the'MSTRG.XXXX' gene id is used.
* In some cases, a novel transcript hasn't been provided a 'reference gene id' yet.
* Some genes have been “merged” due to genomic overlap of transcripts from the different genes.

IsoformSwitchAnalyzeR 's rescue algorithm ensures that 'MSTRG.XXXX' ids are replaced with 'reference gene ids' wherever annotation for the same is available. The result is that only novel genes, those are StringTie identified transcripts that do not overlap with any annotated genes, will still have 'MSTRG.XXXX' ids. 

This helps save data from tens to thousands of genes and can be used for further downstream analysis.

Please install IsoformSwitchAnalyzeR in RStudio befoe you proceed to the commands.

```r
BiocManager::install(IsoformSwitchAnalyzeR)
```

**COMMANDS**

We begin by loading the necessary libraries, reading in the StringTie data, creating a design matrix and combining this into a switchAnalyzeRlist object. The input data will be the 'ballgown' folder that we created while using StringTie. 

It has exon-level, intron-level as well as transcript-level expression measurements for each sample saved in a separate folder.

The switchAnalyzeRlist object is created to specifically contain and summarize all relevant information about the isoforms involved in isoform switches.

Further the gene count matrix is extracted from the switchAnalyzeRlist object.

*Set the working directory to dge_analysis. Also remember to change the paths to the files in the commands below so as to suit your requirements.*

```r
# Loading libraries
library(IsoformSwitchAnalyzeR)

# Import StringTie Expression data
stringTieQuant <- importIsoformExpression(
  parentDir = "rna_seq_dge_analysis/quantification",readLength = 50,
  addIsofomIdAsColumn = FALSE)

# Make design matrix
myDesign <- data.frame(
  sampleID = colnames(stringTieQuant$abundance),
  condition = c("Uninfected","SCovR","SCov",
                "Uninfected","SCovR","SCov"))

# Create switchAnalyzeRlist
switchAnalyzeRlist <- importRdata(
  isoformCountMatrix   = stringTieQuant$counts,
  isoformRepExpression = stringTieQuant$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = "rna_seq_dge_analysis/assembly/merged_transcripts/stringtie_merged.gtf",
  fixStringTieAnnotationProblem = TRUE)

# Extract gene count matrix
geneCountMatrix <- extractGeneExpression(
  switchAnalyzeRlist,
  extractCounts = TRUE)

# Getting Gene count matrix as a csv file
write.csv(geneCountMatrix,"rna_seq_dge_analysis/dge_analysis/gene_count_matrix_new.csv")
```

The 'gene_count_matrix_new.csv' will be used for differential gene analysis with DeSeq2. 

-----------------------------------------------------------------

### 6. Differential Gene Expression Analysis <a name="diff_gene"></a>

Now, that we have our gene_counts_matrix_new.csv file ,obtained from IsoformSwitchAnalyzeR, we can proceed with Differential Gene Expression Analysis. The file contains 33000+ entries for gene expression.
There are several packages that can be used for Differential Gene Expression Analysis like Ballgown, Limma-voom, edgeR but we will be using DESeq2.

#### What is DESeq2?

Differential gene expression analysis helps us compare the differences in gene expression between two conditions and these differences can be considered responsible for the biological phenomena we are observing.
For example, when comparing Diseased vs Normal samples, genes that are expressed in higher/lower quantities in Disease group when compared to the Normal group, could be involved in the disease pathology.

It is important to quantify these differences and conduct a statistical analysis to ascertain which genes are actually responsible for systematic changes between the conditions and to eliminate those that are responsible
for within condition variability. DESeq2 is an R package that tests for differential expression by using negative binomial generalized linear models. 

*What is negative binomial distribution?*

A binomial experiment is one which has a fixed number of independent trials with only two outcomes: success and familar. The probaility of success is constant and a random variable Y indicates the number of successes.

A negative binomial experiment is almost the same except that the number of trials is not fixed and the random variable Y is the number of trials needed to make r successes.

At the moment, we are working with count data - the number of sequence fragments that are assigned to each gene for each sample. Count data comprises of integers (positive or zero) and the variance of the counts increases with the mean. 

The count values tend to aggregate towards a small bunch of the range leading to a positive skew distribution or long right tail. Count data is modeled with either a Poisson or a negative binomial distribution because:
* It has zero or positive integers.
* The variance is a function of mean.
Both these attributes match the properties of our count data.

Most biological count data are not well approximated by a Poisson distribution because the variance is either less than the mean, an example of underdispersion, or greater than the mean, an example of overdispersion. The negative binomial distribution is a useful distribution for count data with overdispersion.

*What are Generalized Linear Models?* 

 Well if have a response variable Y and predictor variable x, then the linear model would be something similar to our straight line equation (y = mx +c)

<p align="center">
     y<sub>i</sub>=β<sub>0</sub>+β<sub>1</sub>x<sub>i</sub>+ε<sub>i</sub>
     ε∼N(0,σ)
 </p>

Here we have a systematic part β<sub>0</sub>+β<sub>1</sub>x<sub>i</sub> and *random error ε<sub>i</sub>* which is a random draw from a normal distribution with mean zero and variance σ<sup>2</sup>.

But in case of generalized linear models, this equation is more helpful:

<p align="center">
     μ<sub>i</sub>=β<sub>0</sub>+β<sub>1</sub>x<sub>i</sub>

μ = E(Y|X)
y<sub>i</sub>∼N(μ<sub>i</sub>,σ)
</p>

Here our *response* is a random draw from a normal distribution with mean mu and variance σ2, and not the error term. 

A generalized linear model would look like this:
<p align="center">
     g(μ<sub>i</sub>) = η<sub>i</sub> = α + β<sub>1</sub>X<sub>i1</sub> + β<sub>2</sub>X<sub>i2</sub> + · · · + β<sub>k</sub>X<sub>ik</sub>
     μ<sub>i</sub> ≡ E(Y<sub>i</sub>)
     </p>

GLM has three components:
* Random part: The response variable Y
* Deterministic part: The predictor/explanatory variables Xi1 to Xik
* Link function: g(μi) = ηi

The natural link function for the Negative Binomial is the “log link”, η=log(μ).

#### How does DESeq2 work?
DESeq2 uses the negative binomial distribution with a slightly more stringent approach compared to other methods but maintaining good balance between sensitivity and specificity (reducing both false positives and false negatives).
DESeq2 performs differential expression analysis in the following manner:
* It calculates normalization factors for sequencing depth adjustment. This accounts for differences in library size.
* The dispersion parameters (α) in negative binomial distribution are estimated. 
* It fits GLM (Generalized Linear Model) for each gene.
* The Wald tests or the likelihood ratio tests are performed to identify DE genes.

**COMMANDS**

 We will now use our raw counts present in gene_counts_matrix_new.csv file in DESeq2. Please consider using RStudio for implementing this code or writing any R scripts.

##### 6a.Installing and Loading Libraries

```r

BiocManger::install("DESeq2"); library(DESeq2)
BiocManager::install("IHW"); library(IHW)
BiocManager::install("org.Hs.eg.db"); library(org.Hs.eg.db)
install.packages("pheatmap"); library(pheatmap)
install.packages("RColorBrewer"); library(RColorBrewer)
BiocManger::install("ggplot2"); library(ggplot2)
BiocManger::install("EnhancedVolcano");library(EnhancedVolcano)
BiocManger::install("enrichR");library(enrichR)
BiocManger::install("gprofiler2");library(gprofiler2)

```

##### 6b. Reading the count matrix
Next, we will read our raw count data from the file gene_counts_matrix_new.csv. Retain only those entries that have relevant IDs in the current NCBI annotation, which is provided by the org.HS.eg.db package.
Check for any duplicate entries and assign ENSEMBL IDs as rownames to our filtered data. Create the count data matrix.

Next, read the phenodata. Phenodata is a file that stores relevant characteristics about your samples and experiment such as sample id, replicates, treatment group etc.
Create a new dataframe ("col_data") that has necessary information extracted from the phenodata file but are stored as factors.

Here the treatment groups are: 
* Uninfected: "Uninfected" 
* Infected with SARS-Cov-2 : "SCov2"
* Infected with SARS-Cov-2 and treated with Remedesvir : "SCov2_R"

The phenodata.csv file is available in **the tutorial_data folder in the repository.** Save it in the dge_analysis folder on your system.

```r

#read counts from gene_count_matrix
targets <- read.csv("rna_seq_dge_analysis/dge_analysis/gene_count_matrix_new.csv", header=TRUE)
head(targets)

#Keeping counts only pertaining to relevant ENSEMBL ids
idfound <- targets$gene_id %in% mappedRkeys(org.Hs.egENSEMBL)
targets_n <- targets[idfound,]

#checking duplicated records
isTRUE(duplicated(targets_n$gene_id))
rownames(targets_n) = targets_n$gene_id

#creating counts matrix
countdata <- as.matrix(round(targets_n[,3:8]))
head(countdata, 3)

#reading phenodata
pheno_data = read.csv("rna_seq_dge_analysis/dge_analysis/phenodata.csv")
head(pheno_data)

col_data <- data.frame("Sex" = pheno_data$Sex,"Replicate" = pheno_data$Replicate)
rownames(col_data) <- pheno_data$Sample.Id
col_data$Sex <- as.factor(col_data$Sex)
col_data$Replicate <- as.factor(col_data$Replicate)
col_data$Group <- c("Uninfected","SCov2_R","SCov2","Uninfected","SCov2_R","SCov2")
col_data$Group <- as.factor(col_data$Group)
col_data$Group
head(col_data)

````

##### 6c. Create the DESeqDataset Object
DESeqDataset object is the object class used by the DESeq2 package to store the read counts and the intermediate estimated quantities during statistical analysis. It accepts the gene counts as a matrix, the phenodata as 'colData' and also the design formula.

The design formula expresses the variables which will be used in modeling. It is used to estimate the dispersions and to estimate the log2 fold changes of the model. 

In our case, the two variables are 'Group' indicating the treatment/infection status ("Uninfected","SCov2_R","SCov2") and 'Sex' indicating the gender of the donor of the samples.


```r

#creating DESeqDataSet object
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = col_data,
                                 design = ~ Sex + Group)

```

#### 6d. Filtering low-count genes
We pre-filter low count genes as they wouldn't be providing any useful information for further statistical analysis. Since our number of samples and replicates are small,
we perform a minimal pre-filtering to keep only rows that have at least at least 10 reads total. 

Starting with 21885 genes, after filtering we are left with 18457 genes.

```r

#Filtering
nrow(ddsMat) 
keep <- rowSums(counts(ddsMat)) > 1
table(keep)
dds <- ddsMat[keep,]
nrow(dds) 

```

##### 6e. Performing Differential Gene Expression Analysis
Finally we perform the DE Analysis on our filtered DESeqDataset object. The estimation steps performed by DEseq function. It estimates size factors, dispersion and fits the model and conducts hypothesis testing.

```r

#DE Analysis
dds <- DESeq(dds)

```

#### 6f. Visualization to assess data quality
Before, we move to our results, lets examine our data via plots. To test for differential expression, we operate on raw counts and use discrete distributions. 

However for other downstream analyses such as for visualization or clustering ,it might be useful to work with transformed versions of the count data.

The obvious choice of transformation is the logarithm. Since count values for a gene can be zero in some conditions (and non-zero in others), it is suggested to use pseudocounts, i.e. transformations of the form:
<p align = "center">
     y=log<sub>2</sub>(n+n<sub>0</sub>)
</p>
     where n represents the count values and n<sub>0</sub> is a positive constant.

One alternative provided in DESeq2 package is variance stabilizing transformations (VST). VST serves a rational way of choosing parameters equivalent to n0 above. 

It removes the dependence of the variance on the mean, particularly the high variance of the logarithm of count data when the mean is low.

We begin by transforming the counts using vst function. Next we assess the count matrix, by constructing a heatmap of normalized counts where samples are differentiated based on Replicates and Treatment Groups.

Then we examine sample clustering by calculating sample to sample distances based on variance stabilized data and plot a heat map where samples are differentiated based on Sex and Treatment Groups.

Related to the distance matrix is the PCA plot, which shows the samples in the 2D plane spanned by their first two principal components. This type of plot is useful for visualizing the overall effect of experimental covariates.

```r

#Transforming values
vsd <- vst(dds)
head(assay(vsd), 3)

#Sample Heatmap: Group & Replicate based on Normalized counts
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Group","Replicate")])
pheatmap(assay(ntd)[select,], show_rownames=TRUE, annotation_col=df)


#Sample Heatmap: Group & Sex based on variance stabilized data
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Group, vsd$Sex, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#PCA Plot: Group & Sex based on variance stabilized data
pcaData <- plotPCA(vsd, intgroup=c("Group", "Sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Group, shape=Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

```

The above-mentioned code should produce plots like these:
<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/sample_heatmap1.png" width="800" height=400 alt="Sample heatmap image 1"/>
</p>

<p align="center">
     <b>A heatmap of normalized counts where samples are differentiated based on Replicates and Treatment Groups. </b>
</p>

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/sample_heatmap2.png" width="800" height=400 alt="Sample heatmap image 2"/>
</p>

<p align="center">
     <b>A heat map of variance stabilized data where samples are differentiated based on Sex and Treatment Groups. </b>
</p>

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/sample_pca.png" width="800" height=400 alt="Principal Component Analysis"/>
</p>

<p align="center">
     <b>Principal Component Analysis </b>
</p>

##### 6g. Extracting Results from DESeqDataset

Now returning to 'dds', our DESeqDataset object that contains our information about DE genes, we will extract results from it that suit our analysis. 

Our design formula ~Group + Sex indicates that the expression of genes in our samples is affected by two covariates: Group (Treatment/Infection Status) and Sex (Gender of the donor).

We are interested in the changes in the gene expression observed in the 'Group' factor. Our covariate 'Group' consists of three categories: "Uninfected", "SCov2" , "SCov2_R".
We want two comparisons to assess DE genes:
* set 1: SCov2 vs Uninfected OR "Samples infected with SARS-Cov-2" vs "Uninfected"
* set 2: SCov2 vs SCov2_R  OR "Samples infected with SARS-Cov-2" vs "Samples infected with SARS-Cov-2 and treated with Remedesvir"

*Please note what set 1 and set 2 refers to because the code from here on will contain variables with these labels. Eg: res1 means results for set 1 : SCov2 vs Uninfected  while res2 means results for set 2 : SCov2 vs SCov2_R*

We start with extracting results for set 1: SCov2 vs Uninfected. We mention the contrasts as contrast = c("Group","SCov2","Uninfected"). A contrast is a linear combination of estimated log2 fold changes, which can be used to test if differences between groups are equal to zero.

To generate more accurate log2 foldchange estimates, DESeq2 allows for the shrinkage of the LFC estimates toward zero when the information for a gene is low i.e. the gene has low counts and/or high variance. 

LFC shrinkage is done by lfcShrink function and we use the shrinkage package 'ashr' which implements Adaptive Shrinkage using using Empirical Bayes.

To optimize the power of our statistical analysis, we use Independent Hypothesis Weighting. It is a is a multiple testing procedure that increases power compared to the method of Benjamini and Hochberg (FDR) by assigning data-driven weights to each hypothesis. 
     
Finally we examine the number of statistically significant DE genes in our results (p-value < 0.1; same as the p-value used by contributors of the dataset). We get 154 DE genes in set 1.

**Results for set 1: SCov2 vs Uninfected**

```r

#Extracting results for set 1: SCov2 vs Uninfected
res1 <- results(dds,contrast = c("Group","SCov2","Uninfected"))
res1
summary(res1)

#Applying Log Fold Shrinkage
res1_shrunken <- lfcShrink(dds, contrast=c("Group","SCov2","Uninfected"), type="ashr", res=res1)
res1_shrunken

#Using Independent Hypothesis Weigting to get accurate adjusted p-values
ihWRes1 = ihw(res1_shrunken,alpha=0.1)
summary(ihWRes1)
table(ihWRes1$padj <= 0.1) #154 DE genes

```

Similarly, we extract results for set 2: "SCov2 vs SCov2_R". We set the contrasts as contrast = c("Group","SCov2","SCov2_R"). We have 131 DE genes in set 2.

**Results for set 2: "SCov2 vs SCov2_R"**

```r

#Extracting results for set2: SCov2 vs SCov2_R 
res2 <- results(dds,contrast = c("Group","SCov2","SCov2_R"))
res2
summary(res2)

#Applying Log Fold Shrinkage
res2_shrunken <- lfcShrink(dds, contrast=c("Group","SCov2","SCovR"), type="ashr", res=res2)
res2_shrunken

#Using Independent Hypothesis Weigting to get accurate adjusted p-values
ihWRes2 = ihw(res2_shrunken,alpha=0.1)
summary(ihWRes2)
table(ihWRes2$padj <= 0.1) #131 DE gene

```

##### 6h. Annotation

Let's provide some basic annotation to our DE genes in set 1 and set 2 (adjusted p-value < 0.1). Annotation is nothing but features or characteristics that help us identify the genes better such as database IDs, gene symbols, chromosome location, gene type etc. 
     
Here we are adding two features RefSeq IDs and Gene Symbols.

```r

#getting RefSeqIDs
egENSEMBL2EG <- toTable(org.Hs.egENSEMBL2EG)
head(egENSEMBL2EG)

#RefSeqIDs for DE genes in set 1
m1 <- match(rownames(ihWRes1), egENSEMBL2EG$ensembl_id)
ihWRes1$entrez <- egENSEMBL2EG$gene_id[m1]

#RefSeqIDs for DE genes in set 2
m2 <- match(rownames(ihWRes2), egENSEMBL2EG$ensembl_id)
ihWRes2$entrez <- egENSEMBL2EG$gene_id[m2]

#adding gene symbol annotation
egSYM <- toTable(org.Hs.egSYMBOL)
head(egSYM)

#Gene symbols for DE genes in set 1
m3 <- match(ihWRes1$entrez, egSYM$gene_id)
ihWRes1$symbol = egSYM$symbol[m3]

#Gene symbols for DE genes in set 2
m4 <- match(ihWRes1$entrez, egSYM$gene_id)
ihWRes2$symbol = egSYM$symbol[m4]

#remove unnecessary variables
rm(egENSEMBL2EG,egSYM,m1,m2,m3,m4)

```

##### 6i. Visualizing Results with Volcano Plots

Lets create volcanoplot for our DE genes in set 1 (SCov2 vs Uninfected, Fig 4) and set 2 (SCov2 vs SCov2_R, Fig 5). A volcano plot is a type of scatterplot that shows statistical significance (P value) versus magnitude of change (fold change). 

It helps identify genes with large fold changes that are also statistically significant. And these may be the most biologically significant genes.  

Since, we are dealing with small number of samples and replicates we may not see many statistically significant genes here. 

```r

#Volcano Plot for set1 
EnhancedVolcano(ihWRes1,
                lab = ihWRes1$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'SARS-CoV-2 versus uninfected',
                pCutoff = 0.1,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)

#Volcano plot for set 2
EnhancedVolcano(ihWRes2,
                lab = ihWRes2$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'SARS-CoV-2 versus SARS-CoV-2 + remdesivir',
                pCutoff = 0.1,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)

```

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/set1_volcanoplot.png" width="800" height=400 alt="Volcano plot image 1"/>
</p>

<p align="center">
     <b>Volcano Plot of Differentially Expressed Genes in SARS-Cov-2 infected samples vs Uninfected Samples. </b>
</p>

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/set2_volcanoplot.png" width="800" height=400 alt="Volcano plot image 2"/>
</p>

<p align="center">
     <b>Volcano Plot of Differentially Expressed Genes in SARS-Cov-2 infected samples vs SARS-Cov-2 infected samples treated with Remdesivir </b>
</p>

##### 6j. Extracting Statistically significant Differentially Expressed Genes

The authors/contributors of this dataset and related manuscript have described statistical significance as adjusted p-value < 0.1 & LogFoldChange LFC > |1|. We will use this criteria to extract our DE genes.

With this, we will have two dataframes (set 1 and 2):
* SCov2vsUninfected_sig : 105 statistically significant DE genes
* SCov2vsSCov2_R_sig: 75 statistically significant DE genes

**Statistically significant DE genes for set 1**

```r

#Extracting genes considered significant (adjusted p-value < 0.1) in set 1
resSig_1 <- subset(ihWRes1, padj < 0.1)
dim(resSig_1) #154 8

#Top down regulated genes in set 1
head(resSig_1[ order(resSig_1$log2FoldChange), ])
#Top up-regulated genes in set 1
head(resSig_1[ order(resSig_1$log2FoldChange, decreasing = TRUE), ])

#Extracting genes with adjusted p-value < 0.1 and LFC > |1| in set 1
SCov2vsUninfected_sig = subset(resSig_1, log2FoldChange > 1 | log2FoldChange < -1 )
dim(SCov2vsUninfected_sig) #105 8

```

**Statistically Significant DE genes for set 2**

```r

#Extracting genes considered significant (adjusted p-value < 0.1) in set 2
resSig_2 <- subset(ihWRes2, padj < 0.1)
dim(resSig_2) #131 8

#Top down-regulated genes in set 2
head(resSig_2[ order(resSig_2$log2FoldChange), ])
#Top up-regulated genes in set 2
head(resSig_2[ order(resSig_2$log2FoldChange, decreasing = TRUE), ])

#Extracting genes with adjusted p-value < 0.1 and LFC > |1| in set 2
SCov2vsSCov2_R_sig = subset(resSig_2, log2FoldChange > 1 | log2FoldChange < -1 )
dim(SCov2vsSCov2_R_sig) #75 8

```

##### 6j. Pathway Enrichment with enrichR
Pathway enrichment analysis helps us utilise the current knowledge of genes and biological processes and compare our DE gene list with them. It identifies biological pathways that are enriched in a gene list more than would be expected by chance.

We will first use enrichR to conduct pathway enrichment analysis specifically for COVID-19-related disease terms. When you load the enrichR package (library(enrichr)) , 
you should see the following messages:

```r

   #Welcome to enrichR
   #Checking connection ... 
   #Enrichr ... Connection is Live!
   #FlyEnrichr ... Connection is available!
   #WormEnrichr ... Connection is available!
   #YeastEnrichr ... Connection is available!
   #FishEnrichr ... Connection is available!

```

enrichR is a web-based tool providing various types of visualization summaries of collective functions of gene lists as well as an alternative approach to rank enriched terms. 

Since, it is a web-based tool, it is important to check for 'Connection is Live' messages so that you can proceed with your analysis. If the server is down, you might get an error code and you will have to wait till the website is live again.

To perform enrichment analysis, we first set our server as 'Enrichr' since we are working on human genes. Next we get gene symbols of our DE genes in both sets: SCov2 vs Uninfected and SCov2 vs SCov2_R.

Set the database to 'COVID-19_Related_Gene_Sets_2021'. We perform separate enrichment analysis for each set using the enrichr() function. Finally we plot our results using plotEnrich() function.

```r

#Gettting DE Gene symbols for set 1 and set 2
set1_names <- SCov2vsUninfected_sig$symbol
set2_names <- SCov2vsSCov2_R_sig$symbol

#Setting the server to Enrichr for human genes
setEnrichrSite("Enrichr")

#Selecting COVID-19 related Database/Gene Set
dbs <- "COVID-19_Related_Gene_Sets_2021"

#Performing pathway enrichment
enriched_1 <- enrichr(set1_names, dbs)
enriched_2 <- enrichr(set2_names, dbs)

#Examining Results
head(enriched_1[["COVID-19_Related_Gene_Sets_2021"]])
head(enriched_2[["COVID-19_Related_Gene_Sets_2021"]])

#Plotting results
plotEnrich(enriched_1[["COVID-19_Related_Gene_Sets_2021"]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")
plotEnrich(enriched_2[["COVID-19_Related_Gene_Sets_2021"]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")

```
The first two plots show enrichment for first 20 terms that are ranked by p-value. 

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/set1_enrichrplot1.png" width="800" height=400 alt="Set1 enrich plot 1"/>
</p>

<p align="center">
     <b>Enrich Plot of Covid-19 Disease Terms in SARS-Cov-2 infected samples vs Uninfected Samples. </b>
</p>

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/set2_enrichrplot1.png" width="800" height=400 alt="Set 2 enrich plot 1"/>
</p>

<p align="center">
     <b>Enrich Plot of Covid-19 Disease Terms in SARS-Cov-2 infected samples vs SARS-Cov-2 infected samples treated with Remdesivir. </b>
</p>

You can extract terms of your choice from the dataframe (eg: enriched_1) and plot the results.
<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/set1_enrichrplot2.png" width="800" height=400 alt="Custom Set 1 enrich plot 2"/>
</p>

<p align="center">
     <b>Custom Enrich Plot of Covid-19 Disease Terms in SARS-Cov-2 infected samples vs Uninfected Samples. </b>
</p>

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/set2_enrichrplot2.png" width="800" height=400 alt="Custom Set 2 enrich plot 2"/>
</p>

<p align="center">
     <b>Custom Enrich Plot of Covid-19 Disease Terms in SARS-Cov-2 infected samples vs SARS-Cov-2 infected samples treated with Remdesivir. </b>
</p>

##### Pathway enrichment with gprofiler2

gprofiler2 provides an R interface to the widely used web toolset g:Profiler (https://biit.cs.ut.ee/gprofiler). It performs functional enrichment analysis and visualization of gene lists, converts gene/protein/SNP identifiers to numerous namespaces, and maps orthologous genes across species. It relies primarily on ENSEMBL databases.

To begin the analysis, we start with extracting up-regulated and downregulated genes and arrange them in decreasing order of Log Fold Change Values. This is done for each set. 

Next, we perform the enrichment analysis using the gost() function. Here we set certain parameters:
* Databases : sources=c("GO","KEGG","REAC","WP")
* Organism: organism = "hsapiens"
* P-value threshold: user_threshold = 0.05
* Multiple testing correction method: correction_method = "fdr"

Our results are stored in a dataframe. To plot our results, we convert the p-values to a negative log10 scale. Finally we plot the first twenty entries in our results dataframe using ggplot2. Note that these are first twenty entries, they are not ranked by p-values.

Lastly, you can extract terms of your choice, that are relevant to your current study, and plot results accordingly.

```r

#SET1
#Extracting up-regulated and down-regulated genes in set 1 and arranging them in decreasing LFC values
set1_up = subset(SCov2vsUninfected_sig, log2FoldChange > 1) 
set1_up_ordered = set1_up[order(set1_up$log2FoldChange, decreasing = TRUE),]

set1_down = subset(SCov2vsUninfected_sig, log2FoldChange < -1) 
set1_down_ordered = set1_down[order(set1_down$log2FoldChange),]

#Performing Pathway enrichment for set 1
set1_gp = gost(list("up-regulated" = set1_up_ordered$symbol,"down-regulated" = set1_down_ordered$symbol),
                       organism = "hsapiens",ordered_query = TRUE,user_threshold = 0.05, correction_method = "fdr", 
                       multi_query = FALSE, evcodes = TRUE, sources=c("GO","KEGG","REAC","WP"))

#Transforming p value on negative log10 scale
set1_gp_df = set1_gp$result
set1_gp_df$log_pval = -1 * log10(set1_gp_df$p_value)
head(set1_gp_df)


#Pathway Enrichment plot for first 20 entries in gProfiler result
#Note the entries are not ordered by value
ggplot(set1_gp_df[0:20,]) +
  # set overall appearance of the plot
  theme_bw() +
  # Define the dependent and independent variables
  aes(x = term_name , y = log_pval, fill = source) +
  # From the defined variables, create a vertical bar chart
  geom_col() +
  scale_fill_brewer(palette = "Spectral") +
  # Set main and axis titles
  ggtitle("Enriched in SARS-Cov-2 vs Uninfected") +
  xlab("Pathways") +
  ylab("Negative Log 10 Adjusted P-value") +
  # Add a line showing the alpha = 0.01 level
  geom_hline(yintercept = -log10(0.001), linetype='dotted',size = 1) +
  geom_hline(yintercept = -log10(0.01), linetype='dotted',size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype='dotted',size = 1) +
  # Flip the x and y axes
  coord_flip()


#SET2
#Extracting up-regulated and down-regulated genes in set 1 and arranging them in decreasing LFC values
set2_up = subset(SCov2vsSCov2_R_sig, log2FoldChange > 1) 
set2_up_ordered = set2_up[order(set2_up$log2FoldChange, decreasing = TRUE),]

set2_down = subset(SCov2vsSCov2_R_sig, log2FoldChange < -1) 
set2_down_ordered = set2_down[order(set2_down$log2FoldChange),]


set2_gp = gost(list("up-regulated" = set2_up_ordered$symbol,"down-regulated" = set2_down_ordered$symbol),
                       organism = "hsapiens",ordered_query = TRUE,user_threshold = 0.05, correction_method = "fdr", 
                       multi_query = FALSE, evcodes = TRUE, sources=c("GO","KEGG","REAC","WP"))

#Transforming p value on negative log10 scale
set2_gp_df = set2_gp$result
set2_gp_df$log_pval = -1 * log10(set2_gp_df$p_value)
head(set2_gp_df)


#Pathway Enrichment plot for first 20 entries in gProfiler result
#Note the entries are not ordered by value
ggplot(set2_gp_df[0:20,]) +
  # set overall appearance of the plot
  theme_bw() +
  # Define the dependent and independent variables
  aes(x = term_name , y = log_pval, fill = source) +
  # From the defined variables, create a vertical bar chart
  geom_col() +
  scale_fill_brewer(palette = "Spectral") +
  # Set main and axis titles
  ggtitle("Enriched in SARS-Cov-2 vs SARS-CoV-2 + remdesivir") +
  xlab("Pathways") +
  ylab("Negative Log 10 Adjusted P-value") +
  # Add a line showing the alpha = 0.01 level
  geom_hline(yintercept = -log10(0.001), linetype='dotted',size = 1) +
  geom_hline(yintercept = -log10(0.01), linetype='dotted',size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype='dotted',size = 1) +
  # Flip the x and y axes
  coord_flip()

```r

The first two plots show enrichment for first 20 terms. 

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/set1_gprofilerplot1.png" width="800" height=400 alt="Set1 gprofiler plot 1"/>
</p>

<p align="center">
     <b>GProfiler Plot for Gene Enrichment Analysis in SARS-Cov-2 infected samples vs Uninfected Samples. </b>
</p>

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/set2_gprofilerplot1.png" width="800" height=400 alt="Set 2 gprofiler plot 1"/>
</p>

<p align="center">
     <b>GProfiler Plot of Gene Enrichment Analysis in SARS-Cov-2 infected samples vs SARS-Cov-2 infected samples treated with Remdesivir. </b>
</p>

You can extract terms of your choice from the dataframe (eg: set1_gp_df) and plot the results.
<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/set1_gprofilerplot2.png" width="800" height=400 alt="Custom Set 1 gprofiler plot 2"/>
</p>

<p align="center">
     <b>Custom GProfiler Plot of Gene Enrichment Analysis in SARS-Cov-2 infected samples vs Uninfected Samples. </b>
</p>

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/images/set2_gprofilerplot2.png" width="800" height=400 alt="Custom Set 2 gprofiler plot 2"/>
</p>

<p align="center">
     <b>Custom GProfiler Plot of Gene Enrichment Analysis in SARS-Cov-2 infected samples vs SARS-Cov-2 infected samples treated with Remdesivir. </b>
</p>

-------------------------------------------------------------------------------
## Citations <a name="citations_list"></a>

**[1]** Müller JA, Groß R, Conzelmann C, et al. "**SARS-CoV-2 infects and replicates in cells of the human endocrine and exocrine pancreas.**" 
        Nat Metab. 2021;3(2):149-165. doi:10.1038/s42255-021-00347-1 PMID: 33536639. [[Research paper](https://pubmed.ncbi.nlm.nih.gov/33536639/)]

**[2]** Andrews S. (2010)."**FastQC: a quality control tool for high throughput sequence data.**" [[Source Code:](http://www.bioinformatics.babraham.ac.uk/projects/fastqc)]
        
**[3]** Bolger, A. M., Lohse, M., & Usadel, B. (2014). "**Trimmomatic: A flexible trimmer for Illumina Sequence Data.**"
        Bioinformatics, btu170. PMID: 24695404 [[Research paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/)]

**[4]** Kim, D., Paggi, J.M., Park, C. et al "**Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype.**"
        Nat Biotechnol 37, 907–915 (2019). https://doi.org/10.1038/s41587-019-0201-4 [[Research paper](https://www.nature.com/articles/s41587-019-0201-4)]

**[5]** Kovaka S, Zimin AV, Pertea GM, Razaghi R, Salzberg SL, Pertea M. "** Transcriptome assembly from long-read RNA-seq alignments with StringTie2**"
        Genome Biology 20, 278 (2019), doi:10.1186/s13059-019-1910-1  [[Research paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1910-1)]

**[6]** Kevin L Howe, Premanand Achuthan, James Allen, Jamie Allen, Jorge Alvarez-Jarreta, M Ridwan Amode, Irina M Armean, Andrey G Azov, Ruth Bennett, Jyothish Bhai, Konstantinos Billis,
        Sanjay Boddu, Mehrnaz Charkhchi, Carla Cummins, Luca Da Rin Fioretto, Claire Davidson, Kamalkumar Dodiya, Bilal El Houdaigui, Reham Fatima, Astrid Gall, Carlos Garcia Giron, 
        Tiago Grego, Cristina Guijarro-Clarke, Leanne Haggerty, Anmol Hemrom, Thibaut Hourlier, Osagie G Izuogu, Thomas Juettemann, Vinay Kaikala, Mike Kay, Ilias Lavidas, Tuan Le, 
        Diana Lemos, Jose Gonzalez Martinez, José Carlos Marugán, Thomas Maurel, Aoife C McMahon, Shamika Mohanan, Benjamin Moore, Matthieu Muffato, Denye N Oheh, Dimitrios Paraschas, 
        Anne Parker, Andrew Parton, Irina Prosovetskaia, Manoj P Sakthivel, Ahamed I Abdul Salam, Bianca M Schmitt, Helen Schuilenburg, Dan Sheppard, Emily Steed, Michal Szpak, 
        Marek Szuba, Kieron Taylor, Anja Thormann, Glen Threadgold, Brandon Walts, Andrea Winterbottom, Marc Chakiachvili, Ameya Chaubal, Nishadi De Silva, Bethany Flint, Adam Frankish, 
        Sarah E Hunt, Garth R IIsley, Nick Langridge, Jane E Loveland, Fergal J Martin, Jonathan M Mudge, Joanella Morales, Emily Perry, Magali Ruffier, John Tate, David Thybert, 
        Stephen J Trevanion, Fiona Cunningham, Andrew D Yates, Daniel R Zerbino, Paul Flicek. "** Ensembl 2021.**"
        Nucleic Acids Res. 2021, vol. 49(1):884–891 PMID: 33137190.  [[Research paper](https://academic.oup.com/nar/article/49/D1/D884/5952199)]

**[7]** Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. "** Twelve years of SAMtools and BCFtools**"
        GigaScience (2021) 10(2) giab008 [33590861] PMID: 33590861  [[Research paper](https://pubmed.ncbi.nlm.nih.gov/33590861/)]

**[8]** Vitting-Seerup K, Sandelin A (2019). "** IsoformSwitchAnalyzeR: Analysis of changes in genome-wide patterns of alternative splicing and its functional consequences.**"
        Bioinformatics. doi: 10.1093/bioinformatics/btz247 [[Research paper](https://academic.oup.com/bioinformatics/article/35/21/4469/5466456)]

**[9]** Love MI, Huber W, Anders S (2014). "** Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.**"
        Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8.  [[Research paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)]

**[10]** Ignatiadis N, Klaus B, Zaugg J, Huber W (2016). "** Data-driven hypothesis weighting increases detection power in genome-scale multiple testing.**"
        Nature Methods. doi: 10.1038/nmeth.3885.   [[Research paper](https://www.nature.com/articles/nmeth.3885)]

**[11]** Carlson M (2019). "** org.Hs.eg.db: Genome wide annotation for Human. R package version 3.8.2. **"

**[12]** Stephens, M. (2016). "** False discovery rates: a new deal.**" Biostatistics, 18:2. 10.1093/biostatistics/kxw041   [[Research paper](https://academic.oup.com/biostatistics/article/18/2/275/2557030)]

**[13]** IXie Z, Bailey A, Kuleshov MV, Clarke DJB., Evangelista JE, Jenkins SL, Lachmann A, Wojciechowicz ML, Kropiwnicki E, Jagodnik KM, Jeon M, & Ma’ayan A. "** Gene set knowledge discovery with Enrichr.**"
        Current Protocols, 1, e90. 2021. doi: 10.1002/cpz1.90 [[Research paper](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.90)]

**[14]** Kolberg L, Raudvere U, Kuzmin I, Vilo J, Peterson H (2020). "** gprofiler2– an R package for gene list functional enrichment analysis and namespace conversion toolset g:Profiler.**"
        F1000Research, 9 (ELIXIR)(709). R package version 0.2.1.  [[Research paper](https://f1000research.com/articles/9-709)]

**[15]** Blighe K, Rana S, Lewis M (2021). "** EnhancedVolcano: Publication-ready volcano plots with enhanced colouring and labeling. R package version 1.12.0.**"
        [[Source Code](https://github.com/kevinblighe/EnhancedVolcano)]

**[16]** Wickham H (2016). "** ggplot2: Elegant Graphics for Data Analysis.**" Springer-Verlag New York. ISBN 978-3-319-24277-4 [[Source Code](https://ggplot2.tidyverse.org.)]

-----------------------------------------------------------------------------------

## License <a name="license_name"></a>

This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/ShrutiBaikerikar/RNASeq_DGE_tutorial/blob/main/LICENSE) file for details


# Copenhagen

PhD course in analyses of genotyping and next-generation sequencing data in medical and population genetics 2015/2016

Tuesday -  Introduction to NGS data:
 - Lecture 4: Allele frequencies, genotypes and SNPs.
 - Computer Exercises IV

Friday – Detecting selection and genetic load:
 - Lecture 9: Detecting natural selection.
 - Computer Exercise IX

## Material

The data has been already downloaded and it is provided in `/ricco/data/matteo/Data`.
These instructions, including all relevant files and scripts, can be found at `/home/matteo/Copenhagen`.
You do not have to do this, but to download and access all the material in this web page locally you should use [git](http://git-scm.com/).
Uncomment it before running it.
```
# git clone https://github.com/mfumagalli/Copenhagen.git
# cd Copenhagen
# git pull # to be sure you have the latest version, if so you should see "Already up-to-date."
```

## Data

As an illustration, we will use 80 BAM files of human samples (of African, European, East Asian, and Native American descent), a reference genome, and putative ancestral sequence.
*To make things more interesting, we have downsampled our data to an average mean depth of 2X!*.

We will also use VCF files for 120 individuals from the same populations.
The human data represents a small genomic region (1MB on chromosome 2) extracted from the 1000 Genomes Project data set.
More information on this project can be found [here](http://www.1000genomes.org/), including their last publication available [here](http://www.nature.com/nature/journal/v526/n7571/full/nature15393.html).

All data is publicly available.
For your curiosity, a pipeline to retrieve such data is provided [here](Files/data.sh).
Again, you do not have to run this, but if you want to download all data locally you need to have 'samtools', 'bgzip' and 'Rscript' installed in your /usr/bin anrun the following command inside the `Copenhagen` folder in the git package (uncomment it before).
```
# bash Files/data.sh
```
`Data` and `Results` folder will be created automatically inside the `Copenhagen` folder.
Data will be saved (but not pushed to git main repository) in `Data` folder.
Additional scripts are be provided in the `Scripts/` folder.

## Preparation

Please set the path for all programs and data we will be using.
```
# mandatory
ANGSD=/home/ricco/github/angsd
NGSTOOLS=/home/ricco/github/ngsTools
NGSDIST=$NGSTOOLS/ngsDist
MS=/home/ricco/bin/msdir/ms
SS=/home/ricco/github/selscan/bin/linux
# optional
NGSADMIX=/home/ricco/github/angsd/misc/NGSadmix
FASTME=/home/ricco/bin/fastme-2.1.5/binaries/fastme-2.1.5-linux64
```
You also need to provide the location of data and sequences:
```
DIR=/home/matteo/Copenhagen
DATA=/ricco/data/matteo/Data
REF=$DATA/ref.fa.gz
ANC=$DATA/anc.fa.gz
```
Finally, create a folder where you will put all the results.
```
mkdir Results
```

## Case study

*MOTIVATION*

Detecting signatures of natural selection in the genome has the twofold meaning of (i) understanding which adaptive processes shaped genetic variation and (ii) identifying putative functional variants.
In case of humans, biological pathways enriched with selection signatures include pigmentation, immune-system regulation and metabolic processes.
The latter may be related to human adaptation to different diet regimes, depending on local food availability (e.g. the case of lactase persistence in dairy-practicing populations).

The human Ectodysplasin A receptor gene, or EDAR, is part of the EDA signalling pathway which specifies prenatally the location, size and shape of ectodermal appendages (such as hair follicles, teeth and glands).
EDAR is a textbook example of positive selection in East Asians (Sabeti et al. 2007 Nature) with tested phenotypic effects (using transgenic mice).

Recently, a genome-wide association study found the same functional variant in EDAR associated to several human facial traits (era shape, chin protusion, ...) in Native American populations (Adhikari et al. Nat Commun 2016).

*HYPOTHESIS*

- Is the functional allele in East Asian at high frequency in other human populations (e.g. Native Americans)?
- Can we identify signatures of natural selection on EDAR in Native Americans?
- Is selection targeting the same functional variant?

*CHALLENGES*
- Admixed population
- Low-depth sequencing data
- ...

*PLAN OF ACTION*

Goal day 1:

- Estimate allele frequencies for tested variant for African, European, East Asian and Native American samples from low-depth sequencing data
Optional:
- Investigate population structure of American samples related to Europeans and Africans
- Select individuals with high Native American ancestry

Goal day 2:

- Perfom a sliding windows scan based on allele frequency differentiation
- Assess statistical significance of selection signatures through simulations
- Test for extended haplotype homozygosity on high-depth sequencing data

Write the paper!

## Agenda

### Tuesday afternoon -  introduction to NGS data:

#### Lecture 1-2.15pm

* Basics of data handling and filtering
* Maximum likelihood and Bayesian estimation
* Genotype likelihoods
* Allele frequencies, SNPs and genotypes calling
* Brief notes on experimental design
* EM algorithm - Imputation - Phasing (Anders)

#### [Practical](Files/day1.md) 2.30-4pm

* Basic data filtering
* Estimation of allele frequencies and SNP calling
* Genotype calling
* Example: estimation of allele frequencies from low-depth sequencing data: the case of EDAR genetic variation in Native Americans

#### Research lecture

* Detecting adaptive evolution in the high-throughput sequencing era

### Friday morning -  detecting selection

#### Lecture 9-10.15pm

* The effect of selection on the genome
* Methods to detect selection signals
* The problem of assessing significance
* Bias introduced by NGS data
* Summary statistics from low-depth data

#### Practical 10:30-12

* Selection scan based on genetic [differentiation](Files/day2a.md) from low-depth data
* Assessing significance through [simulations](Files/day2b.md)
* Selection test based on [haplotype](Files/day2c.md) diversity
* Example: detection of natural selection from low-depth sequencing data and haplotype data: the case of EDAR genetic variation in Native Americans




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

Please make sure you create `Results` folder:
```
mkdir Results
```
Additional scripts are be provided in the `Scripts/` folder.

## Rationale and preparation

The rationale for the examples and exercises in the practical can be found [here](Files/rationale.sh).

Before running the exercises you need to set all correct paths, as shown [here](Files/preparation.sh).

## Agenda

### Tuesday afternoon -  introduction to NGS data:

#### Lecture 1-2.15pm

* Basics of data handling and filtering
* Maximum likelihood and Bayesian estimation
* Genotype likelihoods
* Allele frequencies, SNPs and genotypes calling
* EM algorithm - Imputation - Phasing (Anders)

#### [Practical](Files/day1.md) 2.30-4pm

* Basic data filtering
* Estimation of allele frequencies and SNP calling
* Genotype calling
* Example: identification of allele frequency differentation from low-depth sequencing data: the case of EDAR genetic variation in Native Americans

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




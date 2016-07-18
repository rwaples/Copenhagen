# Copenhagen

PhD course in analyses of genotyping and next-generation sequencing data in medical and population genetics 2015/2016

Tuesday -  introduction to NGS data:
 - Lecture 4: Allele frequencies, genotypes and SNPs
 - Computer Exercises IV (ANGSD, R, realSFS, GATK/SAMtools)

Friday – Detecting selection and genetic load:
 - Lecture 9: Detecting natural selection. 
 - Computer Exercise IX (e.g. ANGSD/realSFS, PBS 1000 G browser)

## Material

The data has been already downloaded and it is provided here `[give link on cluster]`.
To download and access all the material in this web page locally use [git](http://git-scm.com/).
Uncomment it before running it.
```
# git clone https://github.com/mfumagalli/Copenhagen.git
# cd Copenhagen
# git pull # to be sure you have the latest version, if so you should see "Already up-to-date."
```


## Data

As an illustration, we will use 80 BAM files of human samples (of African, European, East Asian, and Native American descent), a reference genome, and putative ancestral sequence.
We will also use VCF files for 180 individuals from the same populations.
The human data represents a small genomic region (1MB on chromosome 11) extracted from the 1000 Genomes Project data set.
More information on this project can be found [here](http://www.1000genomes.org/), including their last publication available [here](http://www.nature.com/nature/journal/v526/n7571/full/nature15393.html).

All data is publicly available.
A pipeline to retrieve such data is provided [here](https://github.com/mfumagalli/Weggis/blob/master/Files/data.sh).
You need to have 'samtools', 'bgzip' and 'Rscript' installed in your /usr/bin to run this.
Uncomment it before running it.
```
# bash Files/data.sh
```
`Data` and `Results` folder will be created automatically inside the `Weggis` folder.
Data will be saved (but not pushed to git main repository) in `Data` folder.
If you do not use the provided script, please make sure you create `Results` folder:
```
mkdir Results
```

Additional scripts are be provided in the `Scripts/` folder.

## Agenda

### Tuesday afternoon -  introduction to NGS data:

#### Lecture 1-2.15pm

* Basics of data handling and [filtering](link.md)
* Maximum likelihood and Bayesian estimation
* Genotype likelihoods
* Allele frequencies, SNPs and genotypes calling
* EM algorithm - Imputation - Phasing (Anders)
* Paper discussion ?

#### Practical 2.30-4pm

* Basic data filtering
* Estimation of allele frequencies and [SNP calling](link_here.md)
* [Genotype calling](link.md)
* Advanced methods to estimate [SFS](link_here.md)
* Exercise: identification of allele frequency differentation from low-depth sequencing data: the case of XXX genetic variation in Native Americans

#### Research lecture

* Population genomics in the high-throughput sequencing era

### Friday morning -  detecting selection

#### Lecture 9-10.15pm

* The effect of selection on the genome
* Methods to detect selection signals
* The problem of assessing significance
* Bias introduced by NGS data
* Summary statistics from low-depth data

#### Practical 10:30-12

* Selection scan based on [genetic differentiation](https://github.com/mfumagalli/Weggis/blob/master/Files/selection.md) from low-depth data
* Assessing significance through [simulations](https://github.com/mfumagalli/Weggis/blob/master/Files/simulation.md)
* Selection test based on [haplotype diversity](https://github.com/mfumagalli/Weggis/blob/master/Files/deep.md) from high-depth data
* Exercise: identification of natural selection from low-depth sequencing data, with admixture assessment and quantification: the case of XXX genetic variation in Native Americans

















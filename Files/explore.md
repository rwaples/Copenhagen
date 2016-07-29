
In this section we will see how to have a general description of the genetic structure our our samples.
This is a preliminary step for all downstream analyses on detecting selection.
We will also take this opportunity to learn some of the programs we are going to use for the rest of this tutorial.

-----------------------

### ANGSD

For most of the examples, we will use the program [ANGSD](http://popgen.dk/wiki/index.php/ANGSD) (Analysis of Next Generation Sequencing Data) developed by Thorfinn Korneliussen and Anders Albrechtsen at the University of Copenhagen.
More information about its rationale and implemented methods can be found [here](http://www.ncbi.nlm.nih.gov/pubmed/25420514).

According to its website *ANGSD is a software for analyzing next generation sequencing data. The software can handle a number of different input types from mapped reads to imputed genotype probabilities. Most methods take genotype uncertainty into account instead of basing the analysis on called genotypes. This is especially useful for low and medium depth data. The software is written in C++ and has been used on large sample sizes. This program is not for manipulating BAM/CRAM files, but solely a tool to perform various kinds of analysis. We recommend the excellent program SAMtools for outputting and modifying bamfiles.*

By the end of this part you will learn:
* how ANGSD handles input and output files
* how to build up a command line
* how to perform SNP and genotype calling
* how to estimate summary statistics and assess population structure taking data uncertainty into account.

### Preparation

Assuming you are on your home folder:
```
cd /gdc_home4/mfuma
```
and you may want to create a new folder
```
mkdir Wednesday
cd Wednesday
```
and create a folder where you will put all your output files:
```
mkdir Results
```
You also need to provide the location of data and reference sequences:
```
DIR=/gdc_home5/groups/bag2016/wednesday
DATA=$DIR/Data
REF=$DATA/hs37d5.fa
ANC=$DATA/hg19ancNoChr.fa
```

Again, we will use 80 BAM files of human samples (of African, European, East Asian, and Native American descent), a reference genome, and putative ancestral sequence.
The human data represents a small genomic region (1MB on chromosome 11) extracted from the 1000 Genomes Project data set, encompassing the FADS gene family.

Also, to make things more interesting, we have downsampled our data to an average mean depth of *2X*.

--------------------------

ANGSD can accept several input files, as described [here](http://popgen.dk/angsd/index.php/Input):

* BAM/CRAM
* Pileup
* Genotype likelihood/probability files
* VCF

First load the correct module
```
module load angsd
```
To see a full list of options in ANGSD type:
```
angsd
```
and you should see something like
```
...
Overview of methods:
        -GL             Estimate genotype likelihoods
        -doCounts       Calculate various counts statistics
        -doAsso         Perform association study
        -doMaf          Estimate allele frequencies
        -doError        Estimate the type specific error rates
        -doAncError     Estimate the errorrate based on perfect fastas
        -doHWE          Est inbreedning per site
        -doGeno         Call genotypes
        -doFasta        Generate a fasta for a BAM file
        -doAbbababa     Perform an ABBA-BABA test
        -sites          Analyse specific sites (can force major/minor)
        -doSaf          Estimate the SFS and/or neutrality tests genotype calling
        -doHetPlas      Estimate hetplasmy by calculating a pooled haploid frequency

        Below are options that can be usefull
        -bam            Options relating to bam reading
        -doMajorMinor   Infer the major/minor using different approaches
        -ref/-anc       Read reference or ancestral genome
        -doSNPstat      Calculate various SNPstat
        many others

For information of specific options type:
        ./angsd METHODNAME eg
                ./angsd -GL
                ./angsd -doMaf
                ./angsd -doAsso etc
                ./angsd sites for information about indexing -sites files
Examples:
        Estimate MAF for bam files in 'list'
                './angsd -bam list -GL 2 -doMaf 2 -out RES -doMajorMinor 1'
```
Be sure you are using this version `-> angsd version: 0.910-65-g86b939f (htslib: 1.3-29-g091c89c) build(Feb 24 2016 09:27:30)` (check the first lines of the output on the screen).

ANGSD can also perform some basic filtering of the data.
These filters are based on:

* quality and depth, see [here](http://www.popgen.dk/angsd/index.php/Filters)
* SNP quality, see [here](http://popgen.dk/angsd/index.php/SnpFilters)
* sites, see [here](http://popgen.dk/angsd/index.php/Sites)

Let us build our command line.
First we need to define input and output files:
```
# angsd -P 4 -b $DIR/ALL_noCHB.bamlist -ref $REF -out Results/ALL \
...
```
with `-b` we give the file including paths to all BAM files we need to analyse.
`-P 4` says the we are going to use 4 threads.
`-ref` specifies the reference sequence.
`-out` states the prefix for all output files that will be generated.

Next we need to define some basic filtering options.
First we define filters based on reads qaulity.
```
# angsd -P 4 -b $DIR/ALL_noCHB.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
...
```
These filters will retain only uniquely mapping reads, not tagged as bad, considering only proper pairs, without trimming, and adjusting for indel/mapping (as in samtools).

Then, we define filters based on mapping and base quality, as well as depth.
```
# angsd -P 4 -b $DIR/ALL_noCHB.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 30 -setMaxDepth 300 -doCounts 1 \
...
```
Specifically here we analyse only reads with a minimum mapping quality of 20, and bases with a minimum quality of 20 (the values are in phred scores).
Also we specify that we are analysing only sites where we have data for half of the individuals (30) and minimum and maximum TOTAL depth of 30 and 300, respectively.
`-doCounts 1` simply forces the calculation of depth.

For more details and examples on how to filtering data with ANGSD see [here](https://github.com/mfumagalli/WoodsHole/blob/master/Files/filtering.md).

Next we are going to illustrate other features in ANGSD by analyses of population structure within our sample.

--------------------

### Population structure

We want to investigate the population structure our samples: PEL,(Peruvians), TSI (Europeans), LWK (Africans), and CHB (East Asians).
In particular we are interested to assess whether our PEL samples are genetically admixed with TSI and LWK.

One solution would be to perform a Principal Component Analysis (PCA), Multidimensional Scaling (MDS) or some clustering based on genetic distances among samples.
Then we can check whether some PEL fall within EUR/LWK or all PEL form a separate clade.

To do this, we first need to assign genotypes (or their associated probabilities).
We now see how to use ANGSD to call genotypes.
The specific option is `-doGeno`:
```
angsd -doGeno
...
-doGeno 0
        1: write major and minor
        2: write the called genotype encoded as -1,0,1,2, -1=not called
        4: write the called genotype directly: eg AA,AC etc
        8: write the posterior probability of all possible genotypes
        16: write the posterior probability of called gentype
        32: write the posterior probability of called gentype as binary
        -> A combination of the above can be choosen by summing the values, EG write 0,1,2 types with majorminor as -doGeno 3
        -postCutoff=0.333333 (Only genotype to missing if below this threshold)
        -geno_minDepth=-1       (-1 indicates no cutof)
        -geno_maxDepth=-1       (-1 indicates no cutof)
        -geno_minMM=-1.000000   (minimum fraction af major-minor bases)
        -minInd=0       (only keep sites if you call genotypes from this number of individuals)

        NB When writing the posterior the -postCutoff is not used
        NB geno_minDepth requires -doCounts
        NB geno_maxDepth requires -doCounts
```

Therefore, if we set `-doGeno 2`, genotypes are coded as 0,1,2, as the number of alternate alleles.
If we want to print the major and minor alleles as well then we set `-doGeno 3`.

To calculate the posterior probability of genotypes we need to define a model.
```
angsd -doPost

...
-doPost 0       (Calculate posterior prob 3xgprob)
        1: Using frequency as prior
        2: Using uniform prior
...
```
`-doPost 1` uses the estimate per-site allele frequency as a prior for genotype proportions, assuming Hardy Weinberg Equilibrium.
When the assumption of HWE is not valid, you can use an estimate of the inbreeding coefficient, for instance calculated using [ngsF](https://github.com/fgvieira/ngsF).

For instance, recalling what previously shown as input/output/filtering options, a command line to call genotypes would be:
```
# angsd -P 4 -b $DIR/ALL_noCHB.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 30 -setMaxDepth 300 -doCounts 1 \
        -doGeno 3 -doPost 1
	...
``` 
This command will not run, but it serves as an illustration of how to set parameters for genotype calling.

For more details and examples on genotype calling with ANGSD see [here](https://github.com/mfumagalli/WoodsHole/blob/master/Files/genocall.md).

--------------------------------

Furthermore, we want to restrict this analysis on a set of putative polymorphic sites (SNPs), as non-variable sites (across all samples) will not carry information regarding population structure or differentiation.

This search can be firstly translated in the estimation of the allele frequency for each site.
In other words, at each site we want to to estimate (or count) how many copies of different alleles (two in case of biallelic SNPs) we observe in our sample (across all sequenced individuals).

ANGSD has an option to estimate **allele frequencies**:

```
angsd -doMaf
...
-doMaf  0 (Calculate persite frequencies '.mafs.gz')
        1: Frequency (fixed major and minor)
        2: Frequency (fixed major unknown minor)
        4: Frequency from genotype probabilities
        8: AlleleCounts based method (known major minor)
        NB. Filedumping is supressed if value is negative
-doPost 0       (Calculate posterior prob 3xgprob)
        1: Using frequency as prior
        2: Using uniform prior
Filters:
        -minMaf         -1.000000       (Remove sites with MAF below)
        -SNP_pval       1.000000        (Remove sites with a pvalue larger)
        -rmTriallelic   0.000000        (Remove sites with a pvalue lower)
Extras:
        -ref    (null)  (Filename for fasta reference)
        -anc    (null)  (Filename for fasta ancestral)
        -eps    0.001000 [Only used for -doMaf &8]
        -beagleProb     0 (Dump beagle style postprobs)
        -indFname       (null) (file containing individual inbreedcoeficients)
NB These frequency estimators requires major/minor -doMajorMinor
```

Therefore, the estimation of allele frequencies requires the specification of how to assign the major and minor alleles (if biallelic).
```
angsd -doMajorMinor
...
        -doMajorMinor   0
        1: Infer major and minor from GL
        2: Infer major and minor from allele counts
        3: use major and minor from a file (requires -sites file.txt)
        4: Use reference allele as major (requires -ref)
        5: Use ancestral allele as major (requires -anc)
        -rmTrans: remove transitions 0
        -skipTriallelic 0
```

Finally, you need to specify which genotype likelihood model to use.
```
angsd -GL
...
        -GL=0:
        1: SAMtools
        2: GATK
        3: SOAPsnp
        4: SYK
        5: phys
        -trim           0               (zero means no trimming)
        -tmpdir         angsd_tmpdir/   (used by SOAPsnp)
        -errors         (null)          (used by SYK)
        -minInd         0               (0 indicates no filtering)

Filedumping:
        -doGlf  0
        1: binary glf (10 log likes)    .glf.gz
        2: beagle likelihood file       .beagle.gz
        3: binary 3 times likelihood    .glf.gz
        4: text version (10 log likes)  .glf.gz
```

We may be interested in looking at allele frequencies only for sites that are actually variable in our sample.
Therefore we want to perform a **SNP calling**.
There are two main ways to call SNPs using ANGSD with these options:
```
        -minMaf         0.000000        (Remove sites with MAF below)
        -SNP_pval       1.000000        (Remove sites with a pvalue larger)
```
Therefore we can consider assigning as SNPs sites whose estimated allele frequency is above a certain threhsold (e.g. the frequency of a singleton) or whose probability of being variable is above a specified value.

For more details and examples on SNP calling with ANGSD see [here](https://github.com/mfumagalli/WoodsHole/blob/master/Files/filtering.md).

--------------------------------

Back to our example, we need to compute genotype posterior probabilities for all samples with ANGSD only on putative polymorphic sites.

First, let us see how to perform a hard SNP/genotype calling in ANGSD, assigning individual genotypes.
```
angsd -P 4 -b $DIR/ALL_noCHB.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 30 -setMaxDepth 300 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-3 \
        -doGeno 3 -doPost 1 &> /dev/null
```
Open the output file:
```
less -S Results/ALL.geno.gz
```
The columns are: chromosome, position, major allele, minor allele, genotypes is 0,1,2 format.

We have also generated estimates of the minor allele frequencies and these are stored in .mafs.gz file:
```
less -S Results/ALL.mafs.gz
```
The columns are: chromosome, position, major allele, minor allele, reference allele, allele frequency, p-value for SNP calling, number of individuals with data.

However, since our data is low-depth, genotypes cannot be assigned with high confidence and therefore we want to use **genotype posterior probabilities** instead, using options `-doGeno 8 -doPost 1`.

Let us build our command line, recalling what we have previously defined.
Print out the called genotypes in allelic format (option 4) and the associated genotype probability (option 16, so the total value is 4+16=20) by using a uniform prior (-doPost 2):
```
angsd -P 4 -b $DIR/ALL_noCHB.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 30 -setMaxDepth 300 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-3 \
        -doGeno 20 -doPost 2 &> /dev/null
```
Open the file:
```
less -S Results/ALL.geno.gz
```
and you will see that now we have probabilities appended to all called genotypes.

As you can see, many of the associated probabilities (<0.90) are low making the assignment of genotypes prone to errors.
Probabilities can be further refined by using a HWE-based prior (option -doPost 1).
Nevertheless any bias in genotype calling will then be downstreamed to all further analyses.

------------------

Specifically, we are now performing a principal component analyses (PCA) without relying on called genotypes, but rather by taking their uncertainty into account.
More specifically, the next program we are going to use (ngsTools) takes as input genotype probabilities in binary format, so we need to specify `-doGeno 32`.
Also, we are using a HWE-based prior with `-doPost 1`.

```
angsd -P 4 -b $DIR/ALL_noCHB.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 30 -setMaxDepth 300 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-3 \
        -doGeno 32 -doPost 1 &> /dev/null
```
Unzip the results but you cannot open it since it is in binary format:
```
gunzip Results/ALL.geno.gz
```

Please unload the ANGSD module now and load the ngsTools and R ones:
```
module unload angsd
module load ngsTools
module load R/3.2.3
```

We are going to use `ngsCovar`, which estimates the covariance matrix between individuals based on genotype probabilities.
Then this matrix will be decomposed into principal componenets which will be investigated for population structure analyses.

If you type
```
ngsCovar
```
you will see a list of possible options.

For instance, we need to define how many sites we have.
To retrieve such value, we can inspect the file with allele frequencies:
```
less -S Results/ALL.mafs.gz
```
How many sites do we have?
```
N_SITES=`zcat Results/ALL.mafs.gz | tail -n+2 | wc -l`
echo $N_SITES
```

Now we can perform a PCA but first estimating the covariance matrix using:
```
ngsCovar -probfile Results/ALL.geno -outfile Results/ALL.covar -nind 60 -nsites $N_SITES -call 0 -norm 0 &> /dev/null
```
with the options `-call 0` meaning that we do not perform genotype calling and `-norm 0` that we are not normalising by allele frequency.
The latter may give more weight to low frequency variants which are harder to estimate.

Look at the output file:
```
less -S Results/ALL.covar
```
which represents a matrix of NxN with N individuals giving the covariance.
Note that this matrix is symmetric.

Finally, we perform an eigenvector decomposition and plot the resulting map:
```
# create a cluster-like file defining the labelling (population) for each sample
Rscript -e 'write.table(cbind(seq(1,60),rep(1,60),c(rep("LWK",20),rep("TSI",20),rep("PEL",20))), row.names=F, sep=" ", col.names=c("FID","IID","CLUSTER"), file="Results/ALL.clst", quote=F)'
# run and plot
Rscript $DIR/Scripts/plotPCA_ngstools.R -i Results/ALL.covar -c 1-2 -a Results/ALL.clst -o Results/ALL.pca.pdf
```
where the parameter `-c 1-2` specifies that we are plotting only the first and second component.
On the screen, you will see a series of numbers.
These are the percentage of explained variance for each component.

Finally, you can open the produced image, by for instance copying it locally to your computer:
```
## run this on your local machine by using your username and password
# scp mfuma@gdcsrv1.ethz.ch:/gdc_home4/mfuma/Wednesday/Results/ALL.pca.pdf .
# open ALL.pca.pdf # or use `evince` on ubuntu
```

Comment on the results.

...thinking...

Therefore, PEL samples appear very close to EUR but separated.
Indeed, among all Latin American populations present in the 1000G project, Peruvians are the least admixed population.
We can either use all these samples or, as described below, compute admixture proportions in order to select a subset of putative Native American (unadmixed) samples.
For the rest of the exercises, we are going to use all PEL sample but optionally you can use only the unadmixed samples.
Note that again the estimation of admixture proportions will be performed taking genotype uncertainty into account

Please unload the ngsTools module and reload the ANGSD one:
```
module unload ngsTools
module load angsd
```

-----------------------

**OPTIONAL**

An alternative approach would be to compute genetic distances first, and then perform a MDS on those.
Here the command lines needed to perform such tasks.

```
# Assuming HWE, without filtering based on probabilities, with SNP calling
angsd -P 4 -b $DIR/ALL_noCHB.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 30 -setMaxDepth 300 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-3 \
        -doGeno 8 -doPost 1 &> /dev/null

#Record how many sites we retrieve.
N_SITES=`zcat Results/ALL.mafs.gz | tail -n+2 | wc -l`
echo $N_SITES

#For plotting purposes, we now create a file with labels indicating the population of interest for each sample.
Rscript -e 'cat(paste(rep(c("LWK","TSI","PEL"),each=20), rep(1:20, 3), sep="_"), sep="\n", file="pops.label")'
cat pops.label

#With [ngsDist](https://github.com/fgvieira/ngsDist) we can compute pairwise genetic distances without relying on individual genotype calls.
module load ngsDist
N_SAMPLES=60
ngsDist -verbose 1 -geno Results/ALL.geno.gz -probs -n_ind $N_SAMPLES -n_sites $N_SITES -labels pops.label -o Results/ALL.dist -n_threads 4
less -S Results/ALL.dist
module unload ngsDist

#From these distances, we can perform a MDS analysis and investigate the population genetic structure of our samples.
tail -n +3 Results/ALL.dist | head -n $N_SAMPLES | Rscript --vanilla --slave Scripts/get_MDS.R --no_header --data_symm -n 4 -m "mds" -o Results/ALL.mds

#We can visualise the pairwise genetic distances in form of a tree (in Newick format).
module load fastme
fastme -D 1 -i Results/ALL.dist -o Results/ALL.tree -m b -n b &> /dev/null
cat Results/ALL.tree
module unload fastme

#We can some R packages to plot the resulting tree.
Rscript -e 'library(ape); library(phangorn); pdf(file="Results/ALL.tree.pdf"); plot(read.tree("Results/ALL.tree"), cex=0.5); dev.off();' &> /dev/null
# again you should copy this .pdf to your local machine in order to visualise it
## run this on your local machine by using your username and password
# scp mfuma@gdcsrv1.ethz.ch:/gdc_home4/mfuma/Wednesday/Results/ALL.tree.pdf .
# open ALL.tree.pdf # or use `evince` on ubuntu
```

------------------------

### ADDITIONAL MATERIAL

### Compute admixture proportions across samples

We use ngsAdmix, which again works on genotype probabilities and not on individual calls.
This is suitable for low-depth data.

ngsAdmix requires genotype likelihoods in BEAGLE format as input.
We can compute these quantities with ANGSD with `-doGlf 2`.
```
angsd -P 4 -b $DIR/ALL_noCHB.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 30 -setMaxDepth 300 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-3 \
        -doGlf 2 &> /dev/null
```
Be sure you are using the latest ANGSD version (>0.900).
Otherwise please unload ngsTools module and reload then ANGSD one.

We assume 3 ancestral populations (Europeans, Africans and Native Americans) making up the genetic diversity of our samples.
Therefore we compute admixture proportions with 3 ancestral components.
```
module load ngsadmix
K=3
ngsdadmix -likes Results/ALL.beagle.gz -K $K -outfiles Results/ALL.admix.K$K -P 4 -minMaf 0.02
module unload ngsadmix
```

Combine samples IDs with admixture proportions and inspect the results.
```
paste ALL_noCHB.bamlist Results/ALL.admix.K$K.qopt > Results/ALL.admix.K$K.txt
less -S Results/ALL.admix.K$K.txt
```

From these quantities we can extract how many samples (and which ones) have a high proportion of Native American ancestry (e.g. >0.90).
We can also plot the individual ancestral proportions for PEL samples.
We need to specify the order of ancestral populations in the admixture file `Results/ALL.admix.K$K.txt`.
In my case these are PEL LWK TSI
```
Rscript Scripts/getUnadmixed.R 0.90 PEL LWK TSI
```
Inspect the results.
```
less -S Results/PEL_unadm.BAMs.txt
# if you want to see the image, again copy it to your local machine
# evince Results/ALL.admix.PEL.pdf
```

Now we have a subset of putative Native American samples.
Then, we could perform all analyses considering only these samples.

------------------------

[HOME](https://github.com/mfumagalli/Weggis)


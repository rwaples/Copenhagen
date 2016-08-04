
In this session you will learn how to do:
* basic post-mapping filtering
* variant calling
* allele frequency estimation
* genotype calling

For most of the examples, we will use the program [ANGSD](http://popgen.dk/wiki/index.php/ANGSD) (Analysis of Next Generation Sequencing Data) developed by Thorfinn Korneliussen and Anders Albrechtsen at the University of Copenhagen.
More information about its rationale and implemented methods can be found [here](http://www.ncbi.nlm.nih.gov/pubmed/25420514).

According to its website *ANGSD is a software for analyzing next generation sequencing data. The software can handle a number of different input types from mapped reads to imputed genotype probabilities. Most methods take genotype uncertainty into account instead of basing the analysis on called genotypes. This is especially useful for low and medium depth data. The software is written in C++ and has been used on large sample sizes. This program is not for manipulating BAM/CRAM files, but solely a tool to perform various kinds of analysis. We recommend the excellent program SAMtools for outputting and modifying bamfiles.*

------------------

#### Preparation

We will use 80 BAM files of human samples (of African, European, East Asians, and Native American descent), a reference genome, and putative ancestral sequence.
The human data represents a small genomic region (3M bp on chromosome 2) extracted from the 1000 Genomes Project data set.

To make things more interesting, we have downsampled our data to an average mean depth of *2X*.

Please set the path for all programs and data we will be using.
As an example these are my paths.
```
SAMTOOLS=/data/data/Software/samtools-1.3
NGSTOOLS=/data/Software/ngsTools
ANGSD=/data/data/Software/angsd
NGSDIST=$NGSTOOLS/ngsDist
NGSADMIX=/data/data/Software/NGSadmix/NGSadmix
FASTME=/data/data/Software/fastme-2.1.4/src/fastme
MS=/data/data/Software/msdir/ms
SS=/data/data/Software/selscan/bin/linux
```
However, these paths have been sym-linked to your /usr/bin so they can be called by simply typing their name, e.g. `angsd`.

You also need to provide the location of data and sequences:
```
DATA=/data/data/tmp/Copenhagen/Data
REF=$DATA/ref.fa.gz
ANC=$ANC/anc.fa.gz
```

--------------------------------------------------------------------------

#### Data filtering

First, we will learn how to build a command line in ANGSD, with the specific example of doing a basic post mapping filtering.

To see a full list of options in ANGSD type:
```
$ANGSD/angsd
```
and you should see something like
```
...
Overview of methods:
	-GL		Estimate genotype likelihoods
	-doCounts	Calculate various counts statistics
	-doAsso		Perform association study
	-doMaf		Estimate allele frequencies
	-doError	Estimate the type specific error rates
	-doAncError	Estimate the errorrate based on perfect fastas
	-HWE_pval		Est inbreedning per site or use as filter
	-doGeno		Call genotypes
	-doFasta	Generate a fasta for a BAM file
	-doAbbababa	Perform an ABBA-BABA test
	-sites		Analyse specific sites (can force major/minor)
	-doSaf		Estimate the SFS and/or neutrality tests genotype calling
	-doHetPlas	Estimate hetplasmy by calculating a pooled haploid frequency

	Below are options that can be usefull
	-bam		Options relating to bam reading
	-doMajorMinor	Infer the major/minor using different approaches
	-ref/-anc	Read reference or ancestral genome
	-doSNPstat	Calculate various SNPstat
	-cigstat	Printout CIGAR stat across readlength
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

We will see later of to perform SNP and genotype calling (and many other things) with ANGSD.

ANGSD can accept several input files, as described [here](http://popgen.dk/angsd/index.php/Input):

* BAM/CRAM
* Pileup
* Genotype likelihood/probability files
* VCF

Here we show how ANGSD can also perform some basic filtering of the data.
These filters are based on:

* quality and depth, see [here](http://www.popgen.dk/angsd/index.php/Filters)
* SNP quality, see [here](http://popgen.dk/angsd/index.php/SnpFilters)
* sites, see [here](http://popgen.dk/angsd/index.php/Sites)

Have a look at our list of BAM files:
```
cat ALL.bamlist
wc -l ALL.bamlist
```

If the input file is in BAM format, the possible options are:
```
$ANGSD/angsd -bam
...
parseArgs_bambi.cpp: bam reader:
	-bam/-b		(null)	(list of BAM/CRAM files)
	-i		(null)	(Single BAM/CRAM file)
	-r		(null)	Supply a single region in commandline (see examples below)
	-rf		(null)	Supply multiple regions in a file (see examples below)
	-remove_bads	1	Discard 'bad' reads, (flag >=256) 
	-uniqueOnly	0	Discards reads that doesnt map uniquely
	-show		0	Mimic 'samtools mpileup' also supply -ref fasta for printing reference column
	-minMapQ	0	Discard reads with mapping quality below
	-minQ		13	Discard bases with base quality below
	-trim		0	Number of based to discard at both ends of the reads
	-trim		0	Number of based to discard at 5 ends of the reads
	-trim		0	Number of based to discard at 3 ends of the reads
	-only_proper_pairs 1	Only use reads where the mate could be mapped
	-C		0	adjust mapQ for excessive mismatches (as SAMtools), supply -ref
	-baq		0	adjust qscores around indels (as SAMtools), supply -ref
	-checkBamHeaders 1	Exit if difference in BAM headers
	-doCheck	1	Keep going even if datafile is not suffixed with .bam/.cram
	-downSample	0.000000	Downsample to the fraction of original data
	-nReads		50	Number of reads to pop from each BAM/CRAMs
	-minChunkSize	250	Minimum size of chunk sent to analyses

Examples for region specification:
		chr:		Use entire chromosome: chr
		chr:start-	Use region from start to end of chr
		chr:-stop	Use region from beginning of chromosome: chr to stop
		chr:start-stop	Use region from start to stop from chromosome: chr
		chr:site	Use single site on chromosome: chr
```

Some basic filtering consists in removing, for instance, read with low quality and/or with multiple hits, and this can be achieved using the parameters ```-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1```.
Let us build such command line.
First we need to define input and output files:
```
# $ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
...
```
with `-b` we give the file including paths to all BAM files we need to analyse.
`-P 4` says the we are going to use 4 threads.
`-ref` specifies the reference sequence.
`-out` states the prefix for all output files that will be generated.

Next we need to define some basic filtering options.
First we define filters based on reads quality.
```
# $ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
#        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
...
```
These filters will retain only uniquely mapping reads, not tagged as bad, considering only proper pairs, without trimming, and adjusting for indel/mapping (as in samtools).
`-C 50` reduces the effect of reads with excessive mismatches, while `-baq 1` computes base alignment quality as explained here ([BAQ](http://samtools.sourceforge.net/mpileup.shtml)) used to rule out false SNPs close to INDELS.

Also, you may want to remove reads with low mapping quality and sites with low quality or covered by few reads (low depth).
Under these circumnstances, the assignment of individual genotypes and SNPs is problematic, and can lead to errors. 

However, it is necessary to know the overall distribution of per-site depth, in order to avoid filtering too many (or few) sites.
We first derive the distribution of quality scores and depth on our data set using ```-doQsDist 1 -doDepth 1```.
```
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL.qc \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 800 -minQ 0 &> /dev/null
```
As an illustration here, ```-maxDepth 800``` corresponds to a per-sample average depth of 10.
This option means that all sites with depth equal or greater than 800 will be binned together.

Have a look at the files generated:
```
ls Results/*
...
-> Output filenames:
		->"Results/ALL.qc.arg"
		->"Results/ALL.qc.qs"
		->"Results/ALL.qc.depthSample"
		->"Results/ALL.qc.depthGlobal"
...
```
```
# counts of quality scores
less -S Results/ALL.qc.qs
# counts of per-sample depth
less -S Results/ALL.qc.depthSample 
wc -l Results/ALL.qc.depthSample # 80
# counts of global depth
less -S Results/ALL.qc.depthGlobal 
```

It is convenient to compute the percentiles of these distributions (and visualize them) in order to make an informative decision on the threshold values we will use for our filtering.
```
Rscript Scripts/plotQC.R Results/ALL.qc 2> /dev/null
```
Have a look at the output files:
```
less -S Results/ALL.qc.info
evince Results/ALL.qc.pdf
``` 

Which values would you choose as sensible thresholds on quality score and global depth (min and max)?
We may also want to remove sites where half of the individual have no data. This is achieved by the -minInd option.
A possible command line would contain the following filtering:
```
...
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 \
...
```
which corresponds to the following scenario:

Parameter | Meaning |
--- | --- |
-minInd 30 | use only sites with data from at least N individuals |
-setMinDepth 210 | minimum total depth |
-setMaxDepth 700 | minimum total depth |

Specifically here we analyse only reads with a minimum mapping quality of 20, and bases with a minimum quality of 20 (the values are in phred scores).
Also we specify that we are analysing only sites where we have data for half of the individuals (30) and minimum and maximum TOTAL depth of 30 and 300, respectively.
`-doCounts 1` simply forces the calculation of depth.

ANGSD can also compute more sophisticated metrics to filter out SNPs, as described [here](http://popgen.dk/angsd/index.php/SnpFilters), mostly based on:
* strand bias
* deviation from HWE
* quality score bias
The strand bias models are described [here](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3532123/). Some examples of strand biases, taken from the previously cited study, can be found [here](https://github.com/mfumagalli/WoodsHole/blob/master/Files/strand_bias.png).
Different models are implemented, as seen [here](https://github.com/mfumagalli/WoodsHole/blob/master/Files/strand_bias_eq.png).

We have now seen how to build a command line in ANGSD with the example of doing a basic post-mapping filtering.

----------------------------

#### Estimation of allele frequencies and SNP calling

We now want to estimate allele frequencies at each site, for the whole sample.
In other words, at each site we want to to estimate (or count) how many copies of different alleles (two in case of biallelic variants) we observe in our sample (across all sequenced individuals).
However with low depth data direct counting of individually assigned genotypes can lead to biased allele frequencies.

ANGSD has an option to estimate **allele frequencies** taking into account data uncertainty from genotype likelihoods:
```
$ANGSD/angsd -doMaf
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
$ANGSD/angsd -doMajorMinor
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
$ANGSD/angsd -GL
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
A description of these different implementation can be found [here](http://www.popgen.dk/angsd/index.php/Genotype_likelihoods).
The GATK model refers to the first GATK paper, SAMtools is somehow more sophisticated (non-independence of errors), SOAPsnp requires a reference sequence for recalibration of quality scores, SYK is error-type specific.
For most applications and data, GATK and SAMtools models should give similar results.

Therefore a possible command line might be:
```
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 40 -setMinDepth 40 -setMaxDepth 400 -doCounts 1 -sites sites.txt \
        -GL 1 -doMajorMinor 4 -doMaf 1 -skipTriallelic 1 &> /dev/null
```
where we specify:
* -GL 1: genotype likelihood model as in SAMtools
* -doMajorMinor 4: force the major allele to be the reference (the minor is inferred)
* -doMaf 1: major and minor are fixed

What are the output files?
```
->"Results/ALL.arg"
->"Results/ALL.mafs.gz"
```
`.args` file is a summary of all options used, while `.mafs.gz` file shows the allele frequencies computed at each site.

Have a look at this file which contains estimates of allele frequency values.
```
zcat Results/ALL.mafs.gz | head
```
and you may see something like
```
chromo	position	major	minor	ref	knownEM	nInd
11	61005992	C	A	C	0.000003	32
11	61005993	C	A	C	0.000003	33
11	61005994	A	C	A	0.000002	33
11	61005995	G	A	G	0.000002	33
11	61005996	C	A	C	0.000002	33
11	61005997	C	A	C	0.000004	32
11	61005998	T	A	T	0.000005	33
11	61005999	G	A	G	0.000005	33
11	61006000	G	A	G	0.000005	34
```


We may be interested in looking at allele frequencies only for sites that are actually variable in our sample.
Therefore we want to perform a **SNP calling**.
There are two main ways to call SNPs using ANGSD with these options:
```
        -minMaf         0.000000        (Remove sites with MAF below)
        -SNP_pval       1.000000        (Remove sites with a pvalue larger)
```
Therefore we can consider assigning as SNPs sites whose estimated allele frequency is above a certain threhsold (e.g. the frequency of a singleton) or whose probability of being variable is above a specified value.

As an illustration, let us call SNPs by computing:
 - genotype likelihoods using GATK method;
 - major and minor alleles inferred from genotype likelihoods;
 - frequency from known major allele but unknown minor;
 - SNPs as those having MAF=>0.01.

```
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -sites sites.txt \
        -GL 2 -doMajorMinor 1 -doMaf 2 -skipTriallelic 1  \
        -minMaf 0.01 &> /dev/null
```

You can have a look at the results:
```
zcat Results/ALL.mafs.gz | head

chromo	position	major	minor	ref	unknownEM	nInd
11	61006040	G	A	G	0.015298	41
11	61006070	G	A	G	0.017361	41
11	61007659	T	G	T	0.014695	31
11	61007783	A	C	A	0.013398	38
11	61007804	T	C	T	0.098128	37
11	61007900	A	C	A	0.011518	40
11	61008088	T	A	T	0.011197	43
11	61008109	C	T	C	0.014019	37
11	61013836	C	A	C	0.020777	40
```

How many SNPs?
```
zcat Results/ALL.mafs.gz | tail -n+2 | wc -l
```

As a general guidance, `-GL 1`, `-doMaf 1/2` and `-doMajorMinor 1` should be the preferred choice when data uncertainty is high.
If interested in analyzing very low frequency SNPs, then `-doMaf 2` should be selected.
When accurate information on reference sequence or outgroup are available, one can use `-doMajorMinor` to 4 or 5.
Also, detecting variable sites based on their probability of being SNPs is generally a better choice than defining a threshold on the allele frequency.
However, various cutoffs and a dedicated filtering should be perform to assess robustenss of your called SNPs.

-------------------------

**EXERCISE**

Try varying the cutoff for SNP calling and record how many sites are predicted to be variable for each scenario.
Identify which sites are not predicted to be variable anymore with a more stringent cutoff (e.g. between a pair of scenario), and plot their allele frequencies.

```
# iterate over some cutoffs
for PV in 0.05 1e-2 1e-4 1e-6
do
        if [ $PV == 0.05 ]; then echo SNP_pval NR_SNPs; fi
        angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL.$PV \
		-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
		-minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 -sites sites.txt \
		-GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
		-SNP_pval $PV &> /dev/null
	echo $PV `zcat Results/ALL.$PV.mafs.gz | tail -n+2 | wc -l`
done
```

A possible output is (your numbers may be different):
```
SNP_pval NR_SNPs
0.05 384
1e-2 344
1e-4 277
1e-6 241
```

Which sites differ from 0.05 and 0.01? What is their frequency?
This script will also print out the first 20 discordant sites (pK.EM is the p-value for the SNP calling test).
```
Rscript -e 'mafs1=read.table(gzfile("Results/ALL.1e-2.mafs.gz"), he=T, strings=F); mafs5=read.table(gzfile("Results/ALL.0.05.mafs.gz"), header=T, stringsAsFact=F); mafs5[!(mafs5[,2] %in% mafs1[,2]),][1:20,]; pdf(file="Results/diff_snpcall.pdf"); par(mfrow=c(1,2)); hist(as.numeric(mafs5[!(mafs5[,2] %in% mafs1[,2]),][,6]), main="Discordant SNPs", xlab="MAF (DAF)"); hist(as.numeric(mafs5[(mafs5[,2] %in% mafs1[,2]),][,6]), main="Concordant SNPs", xlab="MAF"); dev.off();'

evince Results/diff_snpcall.pdf
```

Can you draw some conclusions from these results?
Which frequencies are more difficult to estimate and therefore affect SNP calling?


-----------------------

------------------------

**EXAMPLE**

3) Estimate allele frequencies for SNPs in FADS genes of interest

In ANGSD we can restrict our analyses on a subset of positions of interest using the `-sites` option.
The positions we are looking at are the one found under selection in Inuit, shown [here](https://github.com/mfumagalli/WoodsHole/blob/master/Files/snps_inuit.png): 
- 11 61627960 <br>
- 11 61631510 <br>
- 11 61632310 <br>
- 11 61641717 <br>
- 11 61624414 <br>
- 11 61597212 <br>

The file with these positions need to be formatted as (chromosome positions).
```
> snps.txt
echo 11 61627960 >> snps.txt
echo 11 61631510 >> snps.txt
echo 11 61632310 >> snps.txt
echo 11 61641717 >> snps.txt
echo 11 61624414 >> snps.txt
echo 11 61597212 >> snps.txt
```
Inspect the file.
```
cat snps.txt
```

We need to index this file in order for ANGSD to process it.
```
angsd sites index snps.txt
```

We are interested in calculating the derived allele frequencies, so are using the ancestral sequence to polarise the alleles.
Create new lists of BAM files.
```
head -n 20 ALL.bamlist > LWK.sub.bamlist
tail -n 20 ALL.bamlist > TSI.sub.bamlist
cp Results/PEL_unadm.BAMs.txt PEL.sub.bamlist
```

Run ANGSD to compute allele frequencies.
Here we change the filtering (more relaxed) since we are interested in outputting all sites.
```
for POP in LWK TSI PEL
do
        echo $POP
        angsd -P 4 -b $POP.sub.bamlist -ref $REF -anc $ANC -out Results/$POP \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 1 -setMinDepth 10 -setMaxDepth 500 -doCounts 1 \
                -GL 1 -doMajorMinor 5 -doMaf 1 -skipTriallelic 1 \
                -sites snps.txt &> /dev/null
done
```

Inspect the results.
```
zcat Results/LWK.mafs.gz Results/TSI.mafs.gz Results/PEL.mafs.gz
```

Do you see any allele frequency differentiation?

------------------------------------------

#### Genotype calling

Here we will explore several ways to call genotypes from sequencing data, once SNPs have been assigned.
We need to assign genotypes (or their associated probabilities) to each site for each individual.

The specific option is `-doGeno`:
```
$ANGSD/angsd -doGeno
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
$ANGSD/angsd -doPost

...
-doPost 0       (Calculate posterior prob 3xgprob)
        1: Using frequency as prior
        2: Using uniform prior
...
```
`-doPost 1` uses the estimate per-site allele frequency as a prior for genotype proportions, assuming Hardy Weinberg Equilibrium.
When the assumption of HWE is not valid, you can use an estimate of the inbreeding coefficient, for instance calculated using [ngsF](https://github.com/fgvieira/ngsF).

A typical command for genotype calling assuming HWE is:

```
angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 -sites sites.txt\
	-GL 1 -doMajorMinor 1 -doMaf 2 -skipTriallelic 1 \
	-SNP_pval 1e-3 \
	-doGeno 3 -doPost 1 -postCutoff 0 &> /dev/null
```

Have a look at the output file:
```
less -S Results/ALL.geno.gz
```

How many sites have at least one missing genotype?
```
zcat Results/ALL.geno.gz | grep -1 - | wc -l
```
Why is that?

You can control how to set missing genotype when their confidence is low with `-postCutoff`.
For instance, we can set as missing genotypes when their (highest) genotype posterior probability is below 0.95:

```
angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 -sites sites.txt\
	-GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
	-SNP_pval 1e-2 \
	-doGeno 3 -doPost 1 -postCutoff 0.95 &> /dev/null
```

How many sites do we have in total?
How many sites have at least one missing genotype now?
```
zcat Results/ALL.geno.gz | wc -l
zcat Results/ALL.geno.gz | grep -1 - | wc -l
```

Why are there some many sites with missing genotypes?

The mean depth per sample is around 7/8, therefore genotypes cannot be assigned with very high confidence.

Setting this threshold depends on the mean sequencing depth of your data, as well as your application.
For some analyses you need to work only with high quality genotypes (e.g. measure of proportion of shared SNPs for gene flow estimate), while for others you can be more relaxed (e.g. estimate of overall nucleotide diversity).
We will show later how to accurately estimate summary statistics with low-depth data.

If we use a uniform prior, then the command line requires `-doPost 2`:

```
angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 -sites sites.txt\
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-2 \
        -doGeno 3 -doPost 2 -postCutoff 0.95 &> /dev/null
```

How many sites have at least one missing genotype in this case?
```
zcat Results/ALL.geno.gz | grep -1 - | wc -l
```

Did you expect such difference compared to the case of HWE-based prior?

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

#### Genetic distances

1) Investigate population structure (or clustering) of PEL samples (Peruvians), EUR (Europeans) and LWK (Africans).

One solution would be to perform a PCA, MDS or some clustering based on genetic distances among samples.
Then we can check whether some PEL fall within EUR/LWK or all PEL form a separate clade. 

First, compute genotype posterior probabilities for all samples.
 
```
# Assuming HWE, without filtering based on probabilities, with SNP calling
angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-3 \
        -doGeno 8 -doPost 1 &> /dev/null
```

Record how many sites we retrieve.
```
N_SITES=`zcat Results/ALL.mafs.gz | tail -n+2 | wc -l`
echo $N_SITES
```

Create a file with labels indicating the population of interest for each sample.
```
Rscript -e 'cat(paste(rep(c("LWK","TSI","PEL"),each=20), rep(1:20, 3), sep="_"), sep="\n", file="pops.label")'
cat pops.label
```

With [ngsDist](https://github.com/fgvieira/ngsDist) we can computer pairwise genetic distances without relying on individual genotype calls.
```
ngsDist -verbose 1 -geno Results/ALL.geno.gz -probs -n_ind 60 -n_sites $N_SITES -labels pops.label -o Results/ALL.dist -n_threads 4
less -S Results/ALL.dist
```

We can visualise the pairwise genetic distances in form of a tree.
```
fastme -D 1 -i Results/ALL.dist -o Results/ALL.tree -m b -n b &> /dev/null
cat Results/ALL.tree
```

Plot the tree.
```
Rscript -e 'library(ape); library(phangorn); pdf(file="Results/ALL.tree.pdf"); plot(read.tree("Results/ALL.tree"), cex=0.5); dev.off();' &> /dev/null
evince Results/ALL.tree.pdf
```

Therefore, PEL samples appear very close to EUR.
The next step would be to compute admixture proportions in order to select a subset of putative Native American (unadmixed) samples.



[HOME](https://github.com/mfumagalli/Copenhagen)


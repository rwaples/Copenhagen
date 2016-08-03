
Rationale and preparation

#### Rationale

## CASE STUDY

*MOTIVATION*

Detecting signatures of natural selection in the genome has the twofold meaning of (i) understanding which adaptive processes shaped genetic variation and (ii) identifying putative functional variants.
In case of humans, biological pathways enriched with selection signatures include pigmentation, immune-system regulation and metabolic processes.
The latter may be related to human adaptation to different diet regimes, depending on local food availability (e.g. the case of lactase persistence in dairy-practicing populations).

EDAR... involved in...
EDAR is classic example of positive selection in East Asians with phenotypic effects.
Recently, GWAS study found EDAR associated to... in Native Americans, with the same causal snp.

*HYPOTHESIS*

Is positive selection on EDAR (SNP) targeting Native American populations too?

If so, is selection acting on the same variant as in East Asians?


*CHALLENGES*
- Admixed population
- Low-depth sequencing data
- ...

*PLAN OF ACTION*

- Retrieve genomic data for 1000 Genomes Project for Africans, Europeans, East Asians and Americans (low-depth)
- Investigate population structure of American samples relating to Europeans and Africans [optional]
- Select individuals with high Native American ancestry [optional]
- Perfom a sliding windows scan based on allele frequency differentiation and nucleotide diversity
- Compute allele frequencies for SNPs of interest in Native Americans, Europeans and Africans [optional]
- Assess statistical significance through simulations
- Test for extended haplotype homozygosity based on high-depth sequencing data

------------------------

Objectives:

EDAR


goal day1: estimate allele frequencies of EDAR SNP in all pops (and genotype frequencies)

goal day2: compute PBS, run sims for testing, selscan for iHs?

-----------------------------

#### Preparation


For most of the examples, we will use the program [ANGSD](http://popgen.dk/wiki/index.php/ANGSD) (Analysis of Next Generation Sequencing Data) developed by Thorfinn Korneliussen and Anders Albrechtsen at the University of Copenhagen.
More information about its rationale and implemented methods can be found [here](http://www.ncbi.nlm.nih.gov/pubmed/25420514).

We will use 60 BAM files of human samples (of African, European, and Native American descent), a reference genome, and putative ancestral sequence.
The human data represents a small genomic region (1MB on chromosome 11) extracted from the 1000 Genomes Project data set.

Please set the path for all programs and data we will be using.
As an example these are my paths.
```
SAMTOOLS=/data/data/Software/samtools-1.3
NGSTOOLS=/data/Software/ngsTools
ANGSD=/data/data/Software/angsd
NGSDIST=$NGSTOOLS/ngsDist
NGSADMIX=/data/data/Software/NGSadmix/NGSadmix
FASTME=/data/data/Software/fastme-2.1.4/src/fastme
```
However, these paths have been sym-linked to your /usr/bin so they can be called by simply typing their name, e.g. `angsd`.

Then link all the Data and sequences.
```
DATA=/data/data/tmp/Copenhagen/Data
REF=$DATA/ref.fa.gz
ANC=$ANC/anc.fa.gz
```






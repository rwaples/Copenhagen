
Rationale and preparation

#### Rationale

EDAR


goal day1: estimate allele frequencies of EDAR SNP in all pops (and genotype frequencies)

goal day2: compute PBS, run sims for testing, selscan for iHs?


#### Preparation


For most of the examples, we will use the program [ANGSD](http://popgen.dk/wiki/index.php/ANGSD) (Analysis of Next Generation Sequencing Data) developed by Thorfinn Korneliussen and Anders Albrechtsen at the University of Copenhagen.
More information about its rationale and implemented methods can be found [here](http://www.ncbi.nlm.nih.gov/pubmed/25420514).

We will use 60 BAM files of human samples (of African, European, and Native American descent), a reference genome, and putative ancestral sequence.
The human data represents a small genomic region (1MB on chromosome 11) extracted from the 1000 Genomes Project data set.

## Preparation

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

If you downloaded the data using the provided script, this is what you should specify.
```
REF=Data/ref.fa.gz
ANC=Data/anc.fa.gz
```






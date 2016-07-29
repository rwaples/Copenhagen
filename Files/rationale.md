
## CASE STUDY

*MOTIVATION*

Detecting signatures of natural selection in the genome has the twofold meaning of (i) understanding which adaptive processes shaped genetic variation and (ii) identifying putative functional variants.
In case of humans, biological pathways enriched with selection signatures include pigmentation, immune-system regulation and metabolic processes.
The latter may be related to human adaptation to different diet regimes, depending on local food availability (e.g. the case of lactase persistence in dairy-practicing populations).

Fatty acid desatureses gene family (FADS) include three different genes (1,2,3) located on chromosome 11.
They are important players in regulation of metabolite levels, and are responsible for processing of omega-3 (e.g. from fish oil) and omega-6 (e.g. from land mammals) polyunsaturated fatty acids.
FADS gene have been found to be under positive selection in several human populations, including [Europeans](http://www.ncbi.nlm.nih.gov/pubmed/26595274), [East Asians](http://www.ncbi.nlm.nih.gov/pubmed/26432246), [Africans](http://www.ncbi.nlm.nih.gov/pubmed/22503634), and Greenlandic [Inuit](http://www.ncbi.nlm.nih.gov/pubmed/26383953).
Interestingly, selection seems to act on different genes and variants within FADS in these populations.

Selection acting on FADS genes in Inuit was estimated to start 20k years ago, possibly predating the split between the ancestors of Inuit and Native Americans.
Studies have also reported high allele frequency for some of these FADS variants in Latin American populations.

*HYPOTHESIS*

Is positive selection on FADS genes targeting Native American populations as well?

If so, is selection acting on the same variants as in Inuit?

Can we assess on which haplotypes selection (if any) is acting on?

*CHALLENGES*
- Admixed population
- Low-depth sequencing data
- ...

*PLAN OF ACTION*

- Retrieve genomic data for 1000 Genomes Project for Africans, Europeans, East Asians and Americans (low-depth)
- Investigate population structure of American samples relating to Europeans and Africans
- Select individuals with high Native American ancestry [optional]
- Perfom a sliding windows scan based on allele frequency differentiation and nucleotide diversity
- Compute allele frequencies for SNPs of interest in Native Americans, Europeans and Africans [optional]
- Assess statistical significance through simulations
- Test for extended haplotype homozygosity based on high-depth sequencing data

------------------------

[HOME](https://github.com/mfumagalli/Weggis)




We have some suggestive results pointing towards selection acting in EDAR based on patterns of genetic differentiation against Native Americans and Europeans.
However, we employed low-depth data and our analyses were not based on called genotypes.
Therefore, we did not have the power, for instance, to identify the causal variant (if any) and patterns of haplotype distribution.

It is usual to then perform a targeted deep resequencing of our region of interest, to further refine our selection analysis and highlight putative causal variants.
Here we can even include more samples (as the experimental cost will be lower anyway).

Therefore, we are now using high-depth phased data (in VCF format) for **180** samples from CHB, PEL and another Latin American population (CLM, from Colombia) to assess whether selection signatures are shared across other populations.

--------------------------

As we are using phased data, we are able to perform selection tests based on haplotype diversity and/or homozygosity, as seen during the lecture.
We are going to use the software [selscan](https://github.com/szpiech/selscan), which implements several selection tests (iHS, nSL, XP-EHH and so on).
According to its manual *[these] statistics are designed to use phased genotypes to identify putative regions of recent or ongoing positive selection in genomes. They are all based on the model of a hard selective sweep, where a de novo adaptive mutation arises on a haplotype that quickly sweeps toward fixation, reducing diversity around the locus. If selection is strong enough, this occurs faster than recombination or mutation can act to break up the haplotype, and thus a signal of high haplotype homozygosity can be observed extending from an adaptive locus.*

As an illustration, we are here calculating nSL.
Note that we have already performed a filtering on these VCF files.
First, look at the options in selscan:
```
$SS/selscan --help
```
The basic usage of selscan to compute nSL from VCF file is:
```
$SS/selscan --nsl --vcf Data/PEL.chr2.vcf --out Results/PEL
```
As you can see, unlike for iHS and EHH, we are not required to provide a genetic map since values are computed in windows based on physical distances.

For nSL, the integration is not cut after a certain threshold, but we need to set a value of maximum length for the window (in number of SNPs units) for building the haplotype.
This is set with the option `--max-extend-nsl`.
It is also common to filter out variant with very low frequency.

Therefore our command line might be:
```
$SS/selscan --nsl --vcf Data/PEL.chr2.vcf --out Results/PEL --max-extend-nsl 200 --maf 0.02
```
Have a look at the output file, knowning that the header is:
`< locusID > < physicalPos > < ’1 ’ freq > <sl1 > <sl0 > < unstandardized nSL >`
```
less -S Results/PEL.nsl.out
```

Finally, these unstandardized scores are normalized in allele frequency bins across the entire genome.
This can be achieved by the extra command `norm`:
```
$SS/norm --help
...
--bins <int>: The number of frequency bins in [0,1] for score normalization.
        Default: 100
...
```
Thus, our command would be (note tha ihs and nsl normalisation are equivalent):
```
$SS/norm --ihs --files Results/PEL.nsl.out --bins 20
```
The output file is called `Results/PEL.nsl.out.20bins.norm`.
```
less -S Results/PEL.nsl.out.20bins.norm
```
In theory, we should perform the normalization on genome-wide data (or at least chromosome-wide).

We can calculate nSL for our other populations, CHB and CLM.
```
$SS/selscan --nsl --vcf Data/CHB.chr2.vcf --out Results/CHB --max-extend-nsl 200 --maf 0.02
$SS/norm --ihs --files Results/CHB.nsl.out --bins 20
$SS/selscan --nsl --vcf Data/CLM.chr2.vcf --out Results/CLM --max-extend-nsl 200 --maf 0.02
$SS/norm --ihs --files Results/CLM.nsl.out --bins 20
```
We can plot these results:
```
Rscript Scripts/plotnSL.R Results/PEL.nsl.out.20bins.norm Results/PEL.nsl.pdf
Rscript Scripts/plotnSL.R Results/CHB.nsl.out.20bins.norm Results/CHB.nsl.pdf
Rscript Scripts/plotnSL.R Results/CLM.nsl.out.20bins.norm Results/CLM.nsl.pdf
```
and again copy them to your local machine:
```
# scp mfuma@gdcsrv1.ethz.ch:/gdc_home4/mfuma/Wednesday/Results/PEL.nsl.pdf .
# scp mfuma@gdcsrv1.ethz.ch:/gdc_home4/mfuma/Wednesday/Results/CHB.nsl.pdf .
# scp mfuma@gdcsrv1.ethz.ch:/gdc_home4/mfuma/Wednesday/Results/CLM.nsl.pdf .
# open PEL.nsl.pdf
# open CHB.nsl.pdf
# open CLM.nsl.pdf
```

What is happening here? What is wrong? Why is there not a clear patterns of high/low scores?
Note that in these VCF alleles have been "polarise" according to the reference sequence, and not the ancestral one.

To overcome this issue, either you should use absolute values of nSL and perform a normalisation on the minor allele frequency rather than the non-reference allele frequency.

This exercise stresses the importance of knowing that the reference sequence is just an arbitrary sequence of alleles and should not be considered as any kind of ancestral state.

-----------------------

Another powerful statistic to delect selection from hard sweeps is **XP-EHH**, which measures the differential decay of haplotype homozygosity.
XP-EHH are not standardised based on allele frequency bins and therefore are not sensitive to misspecification of the ancestral state.
We assume that our target population is PEL and our reference population is CHB.

This statistic can be again computed using selscan, provided that we have a genetic map file.
A genetic map for chromosome 11 (based on 1000 Genomes data) is provided here:
```
less -S $DIR/Files/genetic_map_GRCh37_chr11.map
```
However we need to extract only the sites that correspond in our VCF file.
We also need to interpolate over the sites that are not recorded in our genetic map.
A simple R script to do that is here:
```
Rscript $DIR/Scripts/getGenMap.R $DIR/Files/genetic_map_GRCh37_chr11.map $DIR/Data/PEL.chr11.vcf > Results/genetic.map
```

Now we can run XP-EHH giving the resulting file as input.
This may take some time...
```
selscan --xpehh --vcf $DIR/Data/PEL.chr11.vcf --vcf-ref $DIR/Data/CHB.chr11.vcf --map Results/genetic.map --out Results/PEL --threads 4
```
The output file has the header:
`< locusID > < physicalPos > < geneticPos > < popA ’1 ’ freq > < ihhA > < popB ’1 ’ freq > < ihhB > < unstandardized XPEHH >`
```
less -S Results/PEL.xpehh.out
```
Again, we normalise the results (knowing that this should be done genome-wide):
```
norm --xpehh --files Results/PEL.xpehh.out
```
Have a look at the results:
```
less -S Results/PEL.xpehh.out.norm
```
and plot them:
```
Rscript $DIRScripts/plotXPEHH.R Results/PEL.xpehh.out.norm Results/PEL.xpehh.pdf
# scp mfuma@gdcsrv1.ethz.ch:/gdc_home4/mfuma/Wednesday/Results/PEL.xpehh.pdf .
# open PEL.xpehh.pdf
```

Comment on the results.

-----------------------------

**ADDITIONAL MATERIAL**

Finally, we are interested in investigating the haplotype distribution in this region.
Specifically, we want to draw a haplotype network, where all (unique) haplotypes are clustered based on their mutual genetic distance.

First, we convert our VCF files into FASTA files.
```
> Results/FADS.fa
Rscript $DIR/Scripts/vcf2fasta.R $DIR/Data/PEL.fads.vcf PEL Results/PEL.fads.snp >> Results/FADS.fa
Rscript $DIR/Scripts/vcf2fasta.R $DIR/Data/CLM.fads.vcf CLM NULL >> Results/FADS.fa
Rscript $DIR/Scripts/vcf2fasta.R $DIR/Data/CHB.fads.vcf CHB NULL >> Results/FADS.fa
```
Have a look at the resulting file:
```
less -S Results/FADS.fa
```
We are using [pegas](https://bioinformatics.oxfordjournals.org/content/26/3/419.full) package in R to draw haplotype network.
```
Rscript $DIR/Scripts/plotNet.R Results/FADS.fa Results/PEL.fads.snp Results/FADS.pdf 2> /dev/null > Results/FADS.diff
```
Open the plot:
```
# scp mfuma@gdcsrv1.ethz.ch:/gdc_home4/mfuma/Wednesday/Results/FADS.pdf .
# open FADS.pdf
```
Each unique haplotype is represented as a circle whose size is proportional to its frequency.

The variants on the branch separating CHB and PEL common haplotypes are candidates to be causal.
```
less -S Results/FADS.diff
grep -P "V \t XXXVI" Results/FADS.diff | head -n 2 | tail -n 1 > Results/FADS.cause.diff
```

Another useful tool for visualising haplotypes is [PopArt](http://popart.otago.ac.nz/).

------------------------

[HOME](https://github.com/mfumagalli/Copenhagen)





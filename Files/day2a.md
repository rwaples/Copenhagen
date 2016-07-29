
### Allele frequency differentiation

The joint-SFS can be considered as a summary statistics for the level of genetic differentiation between populations.
We have also seen how it can be used as prior information when computing FST (and related metrics) without relying on genotype calling.

Here we see how to compute the 2D-SFS, FST, PBS, and other summary statistics from low-depth data using ANGSD.
Our final goal is to detect signatures of selection in our data.

-------------------------------

One of the most important aspect of data analysis for population genetics is the estimate of the Site Frequency Spectrum (SFS).
SFS records the proportions of sites at different allele frequencies. It can be folded or unfolded, and the latter case implies the use of an outgroup species to define the ancestral state.
SFS is informative on the demography of the population or on selective events (when estimated at a local scale).

We use ANGSD to estimate SFS using on example dataset, using the methods described [here](http://www.ncbi.nlm.nih.gov/pubmed/22911679).
Details on the implementation can be found [here](http://popgen.dk/angsd/index.php/SFS_Estimation).
Briefly, from sequencing data one computes genotype likelihoods (as previously described).
From these quantities ANGSD computes posterior probabilities of Sample Allele Frequency (SAF), for each site.
Finally, an estimate of the SFS is computed.

These steps can be accomplished in ANGSD using `-doSaf 1/2` options and the program `realSFS`.
Again, be sure you are using the latest version of ANGSD (>0.900), by typing `angds` and reading the first lines.
If not, type `module unload ngsTools` and then `module load angsd`

```
angsd -doSaf
...
-doSaf          0
        1: perform multisample GL estimation
        2: use an inbreeding version
        3: calculate genotype probabilities
        4: Assume genotype posteriors as input (still beta)
        -doThetas               0 (calculate thetas)
        -underFlowProtect       0
        -fold                   0 (deprecated)
        -anc                    (null) (ancestral fasta)
        -noTrans                0 (remove transitions)
        -pest                   (null) (prior SFS)
        -isHap                  0 (is haploid beta!)
NB:
          If -pest is supplied in addition to -doSaf then the output will then be posterior probability of the sample allelefrequency for each site
```

The SFS is typically computed for each population separately.
We need to slightly modify the filtering options as now each population has 20 samples.
So now we set `-minInd 20 -setMinDepth 20 -setMaxDepth 200`.
Also, we want to estimate the unfolded SFS and we use a putative ancestral sequence to polarise our alleles (to ancestral and derived states).

We cycle across all populations:
```
for POP in LWK TSI CHB PEL
do
        echo $POP
        angsd -P 4 -b $DIR/$POP.bamlist -ref $REF -anc $ANC -out Results/$POP \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 15 -setMinDepth 20 -setMaxDepth 200 -doCounts 1 \
                -GL 1 -doSaf 1 &> /dev/null
done
```

Have a look at the output file.
```
realSFS print Results/PEL.saf.idx | less -S
```
These values represent the sample allele frequency likelihoods at each site, as seen during the lecture.
So the first value (after the chromosome and position columns) is the likelihood of having 0 copies of the derived allele, the second indicates the probability of having 1 copy and so on.
Note that these values are in log format and scaled so that the maximum is 0.

Can you spot any site which is likely to be variable?
What does this mean? It means that you should look for sites where the highest likelihood does not correspond to allele frequencies of 0 or 1.

The next step would be to use these likelihoods and estimate the overall SFS.
This is achieved by the program `realSFS`.
```
realSFS
        -> ---./realSFS------
        -> EXAMPLES FOR ESTIMATING THE (MULTI) SFS:

        -> Estimate the SFS for entire genome??
        -> ./realSFS afile.saf.idx

        -> 1) Estimate the SFS for entire chromosome 22 ??
        -> ./realSFS afile.saf.idx -r chr22

        -> 2) Estimate the 2d-SFS for entire chromosome 22 ??
        -> ./realSFS afile1.saf.idx  afile2.saf.idx -r chr22

        -> 3) Estimate the SFS for the first 500megabases (this will span multiple chromosomes) ??
        -> ./realSFS afile.saf.idx -nSites 500000000

        -> 4) Estimate the SFS around a gene ??
        -> ./realSFS afile.saf.idx -r chr2:135000000-140000000

        -> Other options [-P nthreads -tole tolerence_for_breaking_EM -maxIter max_nr_iterations -bootstrap number_of_replications]

        -> See realSFS print for possible print options
        -> Use realSFS print_header for printing the header

        ->------------------
        -> NB: Output is now counts of sites instead of log probs!!
        -> NB: You can print data with ./realSFS print afile.saf.idx !!
        -> NB: Higher order SFS can be estimated by simply supplying multiple .saf.idx files!!
        -> NB: Program uses accelerated EM, to use standard EM supply -m 0
```
This command will estimate the SFS for each population:
```
for POP in LWK TSI CHB PEL
do
        echo $POP
        realSFS Results/$POP.saf.idx -P 4 2> /dev/null > Results/$POP.sfs
done
```
The output will be saved in Results/POP.sfs files.

You can now have a look at the output file, for instance for the African (LWK) samples:
```
cat Results/LWK.sfs
```
The first value represent the expected number of sites with derived allele frequency equal to 0, the second column the expected number of sites with frequency equal to 1 and so on.

How many values do you expect?
```
awk -F' ' '{print NF; exit}' Results/LWK.sfs
```
Indeed this represents the unfolded spectrum, so it has 2N+1 values with N diploid individuals.

Why is it so bumpy? Why there are low/zero values?

...thinking...

This maximum likelihood estimation of the SFS should be performed at the whole-genome level to have enough information for the algorithm to converge.
However, for practical reasons, here we could not use large genomic regions.
Also, as we will see later, this region is not really a proxy for neutral evolution so the SFS is not expected to behave neutrally.
Nevertheless, these SFS should be a reasonable prior to be used for estimation of summary statistics.

Let us plot the SFS for each pop using this simple R script.
```
Rscript $DIR/Scripts/plotSFS.R Results/LWK.sfs Results/TSI.sfs Results/CHB.sfs Results/PEL.sfs
## copy the file to your local machine to inspect it
# scp mfuma@gdcsrv1.ethz.ch:/gdc_home4/mfuma/Wednesday/Results/LWK_TSI_CHB_PEL.pdf .
# open LWK_TSI_CHB_PEL.pdf
```

Do they behave like expected?
Which population has more SNPs?
Which population has a higher proportion of common (not rare) variants?
Why is it bumpy?

---------------------------------------

**OPTIONAL**

It is sometimes convenient to generate bootstrapped replicates of the SFS, by sampling with replacements genomic segments.
This could be used for instance to get confidence intervals when using the SFS for demographic inferences.
This can be achieved in ANGSD using:
```
realSFS Results/PEL.saf.idx -bootstrap 10  2> /dev/null > Results/PEL.boots.sfs
cat Results/PEL.boots.sfs
```
This command may take some time.
The output file has one line for each boostrapped replicate.

More examples on how to estimate the SFS with ANGSD can be found [here](https://github.com/mfumagalli/WoodsHole/blob/master/Files/sfs.md).

---------------------------------------

**OPTIONAL**

It is very useful to estimate a **multi-dimensional SFS**, for instance the joint SFS between 2 populations (2D).
This can be used for making inferences on their divergence process (time, migration rate and so on).
However, here we are interested in estimating the 2D-SFS as prior information for our FST/PBS.

An important issue when doing this is to be sure that we are comparing the exactly same corresponding sites between populations.
ANGSD does that automatically and considers only a set of overlapping sites.
The 2D-SFS between all populations and PEL are computed with:
```
for POP in LWK TSI CHB
do
        echo $POP
        realSFS -P 4 Results/$POP.saf.idx Results/PEL.saf.idx 2> /dev/null > Results/$POP.PEL.sfs
done
```

The output file is a flatten matrix, where each value is the count of sites with the corresponding joint frequency ordered as [0,0] [0,1] and so on.
```
less -S Results/LWK.PEL.sfs
```
You can plot it, but you need to define how many samples you have per population.
```
Rscript Scripts/plot2DSFS.R Results/LWK.PEL.sfs 20 20
# copy to your local machine and open it
# scp mfuma@gdcsrv1.ethz.ch:/gdc_home4/mfuma/Wednesday/Results/LWK.PEL.sfs.pdf .
# open Results/LWK.PEL.sfs.pdf
```

You can even estimate SFS with higher order of magnitude.
This command may take some time.
```
realSFS -P 4 Results/LWK.saf.idx Results/TSI.saf.idx Results/PEL.saf.idx 2> /dev/null > Results/LWK.TSI.PEL.sfs
```

------------------------------------
Here we are going to calculate **allele frequency differentiation** using the PBS (population branch statistic) metric.
Again, we can achieve this by avoid genotype calling using ANGSD.
From the sample allele frequencies likelihoods (.saf files) we can estimate PBS using the following pipeline.

Note that here we use the previously calculated SFS as prior information.
Also, PEL is our target population, while CHB and TSI are reference populations.
Therefore, we need to compute also the 2D-SFS between TSI and CHB and PEL (if not already done):
```
for POP in TSI CHB
do
        echo Results/$POP.PEL.sfs
        realSFS -P 4 Results/$POP.saf.idx Results/PEL.saf.idx 2> /dev/null > Results/$POP.PEL.sfs
done
echo Results/TSI.CHB.sfs
realSFS -P 4 Results/TSI.saf.idx Results/CHB.saf.idx 2> /dev/null > Results/TSI.CHB.sfs
```

Look at the output file:
```
less -S Results/CHB.PEL.sfs
```
The output file is a flatten matrix, where each value is the count of sites with the corresponding joint frequency ordered as [0,0] [0,1] and so on.

You can also plot it, but you need to define how many samples you have per population.
```
Rscript $DIR/Scripts/plot2DSFS.R Results/CHB.PEL.sfs 20 20
# copy to your local machine and open it
# scp mfuma@gdcsrv1.ethz.ch:/gdc_home4/mfuma/Wednesday/Results/CHB.PEL.sfs.pdf .
# open CHB.PEL.sfs.pdf
```

-----------------------
The 2D-SFS will be used as prior information for the joint allele frequency probabilities at each site.
From these probabilities we will calculate the population branch statistic (PBS) using the PEL as target population and CHB and TSI as reference populations.
Our goal is to detect selection in PEL in terms of allele frequency differentiation.

Specifically, we are computing a slinding windows scan, with windows of 50kbp and a step of 10kbp.
This can be achieved using the following commands.
1) This command will compute per-site FST indexes
```
realSFS fst index Results/TSI.saf.idx Results/CHB.saf.idx Results/PEL.saf.idx -sfs Results/TSI.CHB.sfs -sfs Results/TSI.PEL.sfs -sfs Results/CHB.PEL.sfs -fstout Results/PEL.pbs &> /dev/null
```
and you can have a look at their values:
```
realSFS fst print Results/PEL.pbs.fst.idx | less -S
```
where columns are: chromosome, position, (a), (a+b) values for the three FST comparisons, where FST is defined as a/(a+b).
Note that FST on multiple SNPs is calculated as sum(a)/sum(a+b).

2) The next command will perform a sliding-window analysis
```
realSFS fst stats2 Results/PEL.pbs.fst.idx -win 50000 -step 10000 > Results/PEL.pbs.txt
```
Have a look at the output file:
```
less -S Results/PEL.pbs.txt
```
The header is:
```
region  chr     midPos  Nsites  Fst01   Fst02   Fst12   PBS0    PBS1    PBS2
```
Where are interested in the column `PB2` which gives the PBS values assuming PEL (coded here as 2) being the target population.
Note that negative PBS and FST values are equivalent to 0.

We are also provided with the individual FST values.
You can see that high values of PBS2 are indeed associated with high values of both Fst02 and Fst12 but not Fst01.

We can plot the results along with the gene annotation.
```
Rscript $DIR/Scripts/plotPBS.R Results/PEL.pbs.txt Results/PEL.pbs.pdf
```
It will also print out the maximum PBS value observed (e.g. `Maximum PBS value: 1.151302`), as this value will be used in the next part.
This script will also plot the PBS variation assuming CHB as the target population, as a comparison.
Again, copy this file to your local machine and inspect it:
```
## run this on your local machine by using your username and password
# scp mfuma@gdcsrv1.ethz.ch:/gdc_home4/mfuma/Wednesday/Results/PEL.pbs.pdf .
# open PEL.pbs.pdf
```
Comment on the results.

...thinking...

There is a local increase in PBS (genetic differentiation in PEL) corresponding to the FADS genes cluster!
Such pattern is not observed assuming CHB being the selected population.
Therefore it looks like selection might target FADS genes in Native Americans as well.
This has not been reported yet.

-------------------------

As an illustration, now we are going to check what we would have got if we had used 'default' values under a standard genotype calling approach.
For the sake of fairness, we are doing a light filtering based on mapping and base quality as well.
Also we assume that the default cutoff for SNP calling is set on a minimum allele frequency of '1&'.
Allele frequencies are simply calculated as proportions of alleles at each site.
Finally, the two alleles are simply defined by allele counts.

Call genotypes under these requirements:
```
angsd -P 4 -b $DIR/ALL.bamlist -ref $REF -out Results/DUMB \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 40 -setMinDepth 10 -setMaxDepth 1000 -doCounts 1 \
        -GL 1 -doMajorMinor 2 -doMaf 8 -skipTriallelic 1 \
        -doGeno 2 -doPost 2 -minMaf 0.01 &>/dev/null
```
and open the results:
```
less -S Results/DUMB.geno.gz
```
as you can see there are many '-1' values which represent unassigned genotypes due to missing data.
We would not be able to perform any analyses with some many missing sites, unless doing some imputation.
We can overcome this by assuming HWE and use an informative prior on genotypes.
This is achieved by using `-doPost 1` option.
Allele frequencies are again simply calculated as proportions of alleles at each site.

Call genotypes using a HWE-based prior:
```
angsd -P 4 -b $DIR/ALL.bamlist -ref $REF -out Results/DUMB \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 40 -setMinDepth 10 -setMaxDepth 1000 -doCounts 1 \
        -GL 1 -doMajorMinor 2 -doMaf 8 -skipTriallelic 1 \
        -doGeno 2 -doPost 1 -minMaf 0.01 &>/dev/null
```
and open the results:
```
less -S Results/DUMB.geno.gz
```

We now calculate PBS value from called genotypes, again in sliding windows in 50kbp with a setp of 10kbp.
You can use this R script to plot the results as well:
```
Rscript $DIR/Scripts/calcPBS.R Results/DUMB.geno.gz Results/DUMB.pbs.pdf
```
and again copy the file to your local machine:
```
# scp mfuma@gdcsrv1.ethz.ch:/gdc_home4/mfuma/Wednesday/Results/DUMB.pbs.pdf .
# open DUMB.pbs.pdf
```

Comment on the results.

...thinking...

As you can see, there is a local increase around the FADS genes but the signal is way reduced compared to the previous analysis.
Note that the scale on the y-axis is the same as before.
Therefore by using a naive approach based on genotype calling we would completely miss the putative selection signal.
For low-depth data, a probabilistic approach that estimates summary statistics taking data uncertainty into account is highly recommeded.

The next step would be to assess whether such pattern of allelic differentiation is expected given a demographic model.

-------------

**OPTIONAL**

We are also interested in assessing whether an increase in allele frequency differentiation is also associated with a change of **nucleotide diversity** in PEL.
Again, we can achieve this using ANGSD by estimating levels of diversity without relying on called genotypes.

The procedure is similar to what done for PBS, and the SFS is again used as a prior to compute allele frequencies probabilities.
From these quantities, expectations of various diversity indexes are compute.
This can be achieved using the following pipeline.

First we compute the allele frequency posterior probabilities and associated statistics (-doThetas) using the SFS as prior information (-pest)
```
for POP in TSI CHB PEL
do
        echo $POP
        angsd -P 4 -b $DIR/$POP.bamlist -ref $REF -anc $ANC -out Results/$POP \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 10 -setMaxDepth 200 -doCounts 1 \
                -GL 1 -doSaf 1 \
                -doThetas 1 -pest Results/$POP.sfs &> /dev/null
done
```
Then we need to index thess file and perform a sliding windows analysis using a window length of 50kbp and a step size of 10kbp.
```
for POP in TSI CHB PEL
do
        echo $POP
        # index files
        thetaStat make_bed Results/$POP.thetas.gz &> /dev/null
        # perform a sliding-window analysis
        thetaStat do_stat Results/$POP.thetas.gz -nChr 1 -win 50000 -step 10000 -outnames Results/$POP.thetas &> /dev/null
done
```
Values in this output file are the sum of the per-site estimates for the whole window.
Output of Tajimas D is currently disabled so we are going to plot the unstandardised difference between Theta Watterson and Pi.
For instance,
```
less -S Results/PEL.thetas.pestPG
```

Finally, plot the sliding windows scan for nuclotide diversity:
```
Rscript $DIR/Scripts/plotSS.R
# copy to your local machine, run from your local machine
# scp mfuma@gdcsrv1.ethz.ch:/gdc_home4/mfuma/Wednesday/Results/PEL.ss.pdf .
# open PEL.ss.pdf
```

Comment on the results.

...thinking...

We do not observe a prominent signal of decay of nucleotide diversity specific to Peruvians.
This may suggest us that the selection signatures might be complex and involve multiple haplotypes.
Therefore, we now need to run simulations to test for significance of such values of PBS.
Next, we will investigate the patterns of haplotype diversity using networks.

---------------------------------------------------------------------------

**OPTIONAL**

We may also be interested in estimating allele frequencies for our SNPs of interest.

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
Run ANGSD to compute allele frequencies.
Here we change the filtering (more relaxed) since we are interested in outputting all sites.
```
for POP in LWK TSI CHB PEL
do
        echo $POP
        angsd -P 4 -b $DIR/$POP.bamlist -ref $REF -anc $ANC -out Results/$POP \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 10 -setMaxDepth 100 -doCounts 1 \
                -GL 1 -doMajorMinor 5 -doMaf 1 -skipTriallelic 1 \
                -sites snps.txt &> /dev/null
done
```

Inspect the results.
```
zcat Results/LWK.mafs.gz Results/TSI.mafs.gz Results/CHB.mafs.gz Results/PEL.mafs.gz
```

Do you see any allele frequency differentiation?


------------------------

[HOME](https://github.com/mfumagalli/Copenhagen)




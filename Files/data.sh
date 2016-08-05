
## Pipeline to download and process the data to be used for this workshop

# set path
# SAMTOOLS=/data/data/Software/samtools-1.3/samtools
SAMTOOLS=samtools
echo Is this your path to samtools? $SAMTOOLS
BGZIP=bgzip
echo Is this your path to bgzip? $BGZIP

VCFLIB=/data/data/Software/vcflib/bin
echo Is this your path to VCFlib? $VCFLIB

# chrom start end
CHROM=2
START=109000000
END=110000000

mkdir Data
mkdir Results

# get unrelated samples iDs
Rscript Scripts/getUnrelated.R

# download BAM files
bash Scripts/getBams.sh $CHROM $START $END
# this creates files and folders in Data/PEL.BAMs/* and TSI and LWK and CHB

# create file with list of BAMs
ls Data/LWK.BAMs/*.bam Data/TSI.BAMs/*.bam Data/CHB.BAMs/*.bam Data/PEL.BAMs/*.bam > ALL.bamlist
ls Data/LWK.BAMs/*.bam Data/TSI.BAMs/*.bam Data/PEL.BAMs/*.bam > ALL_noCHB.bamlist
ls Data/LWK.BAMs/*.bam > LWK.bamlist
ls Data/TSI.BAMs/*.bam > TSI.bamlist
ls Data/CHB.BAMs/*.bam > CHB.bamlist
ls Data/PEL.BAMs/*.bam > PEL.bamlist

# download ancestral sequence
echo Downloading and processing ancestral sequence...
#wget ftp://ftp.ensembl.org/pub/release-65/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e65.tar.bz
#tar xjf homo_sapiens_ancestor_GRCh37_e65.tar.bz
#cp homo_sapiens_ancestor_GRCh37_e65/homo_sapiens_ancestor_$CHROM.fa Data/tmp.fa
#sed "1s/.*/>2/" Data/tmp.fa > Data/anc.fa
#rm Data/tmp.fa
#$BGZIP Data/anc.fa
#$SAMTOOLS faidx Data/anc.fa.gz
#rm -rf homo_sapiens_ancestor_GRCh37_e65*

wget http://popgen.dk/software/download/angsd/hg19ancNoChr.fa.gz
zcat hg19ancNoChr.fa.gz > Data/anc.fa
$BGZIP Data/anc.fa
$SAMTOOLS faidx Data/anc.fa.gz
rm hg19ancNoChr.fa.gz


# download reference sequence
echo Downloading and processing reference sequence...
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz &> /dev/null
zcat human_g1k_v37.fasta.gz > Data/ref.fa 2> /dev/null
$BGZIP Data/ref.fa
$SAMTOOLS faidx Data/ref.fa.gz
rm human_g1k_v37.fasta.gz

# download VCF files
echo Downloading VCF file...
wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr$CHROM.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz &> /dev/null
wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr$CHROM.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi &> /dev/null

echo Processing VCF file...
VCF=ALL.chr$CHROM.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

# only 3Mbp
START=108000000
END=111000000

# whole region for selscan
$VCFLIB/vcffilter -f "VT = SNP" -f "AC > 1" -f "QUAL > 20" -f "DP > 20" -f "DP < 40000" -r $CHROM:$START-$END $VCF > ALL.chr$CHROM.tmp.vcf
grep -v MULTI_ALLELIC ALL.chr$CHROM.tmp.vcf > Data/ALL.chr$CHROM.vcf

# PEL

$VCFLIB/vcfkeepsamples Data/ALL.chr$CHROM.vcf HG01565 HG01566 HG01571 HG01572 HG01577 HG01578 HG01892 HG01893 HG01917 HG01918 HG01920 HG01921 HG01923 HG01924 HG01926 HG01927 HG01932 HG01933 HG01935 HG01938 HG01939 HG01941 HG01942 HG01944 HG01945 HG01947 HG01948 HG01950 HG01951 HG01953 HG01954 HG01961 HG01965 HG01967 HG01968 HG01970 HG01971 HG01973 HG01974 HG01976 HG01977 HG01979 HG01980 HG01982 HG01991 HG01992 HG01995 HG01997 HG02002 HG02003 HG02008 HG02089 HG02090 HG02104 HG02105 HG02146 HG02147 HG02252 HG02253 HG02259 > Data/PEL.chr$CHROM.vcf
sed -i -e 's/\t\./\t0|0/g' Data/PEL.chr$CHROM.vcf

$VCFLIB/vcfkeepsamples Data/ALL.chr$CHROM.vcf NA18524 NA18525 NA18526 NA18527 NA18528 NA18529 NA18530 NA18531 NA18532 NA18533 NA18534 NA18535 NA18536 NA18537 NA18538 NA18539 NA18540 NA18541 NA18542 NA18543 NA18544 NA18545 NA18546 NA18547 NA18548 NA18549 NA18550 NA18552 NA18553 NA18555 NA18557 NA18558 NA18559 NA18560 NA18561 NA18562 NA18563 NA18564 NA18565 NA18566 NA18567 NA18569 NA18570 NA18571 NA18572 NA18573 NA18574 NA18575 NA18576 NA18577 NA18579 NA18580 NA18582 NA18583 NA18591 NA18592 NA18593 NA18595 NA18596 NA18597 > Data/CHB.chr$CHROM.vcf
sed -i -e 's/\t\./\t0|0/g' Data/CHB.chr$CHROM.vcf

$VCFLIB/vcfkeepsamples Data/ALL.chr$CHROM.vcf HG01112 HG01113 HG01119 HG01121 HG01122 HG01124 HG01125 HG01130 HG01131 HG01133 HG01134 HG01136 HG01137 HG01139 HG01140 HG01142 HG01148 HG01149 HG01250 HG01251 HG01253 HG01254 HG01256 HG01257 HG01259 HG01260 HG01269 HG01271 HG01272 HG01274 HG01275 HG01277 HG01278 HG01280 HG01281 HG01284 HG01341 HG01342 HG01344 HG01345 HG01347 HG01348 HG01350 HG01351 HG01353 HG01354 HG01356 HG01357 HG01359 HG01360 HG01362 HG01363 HG01365 HG01366 HG01369 HG01372 HG01374 HG01375 HG01377 HG01378 > Data/CLM.chr$CHROM.vcf
sed -i -e 's/\t\./\t0|0/g' Data/CLM.chr$CHROM.vcf

# EDAR is our target gene
# chr2:109,510,927-109,605,828
START=109510927
END=109605828
GENE=edar

$VCFLIB/vcffilter -f "VT = SNP" -f "AC > 1" -f "QUAL > 20" -f "DP > 20" -f "DP < 40000" -r $CHROM:$START-$END $VCF > ALL.$GENE.tmp.vcf
grep -v MULTI_ALLELIC ALL.$GENE.tmp.vcf > Data/ALL.$GENE.vcf

$VCFLIB/vcfkeepsamples Data/ALL.$GENE.vcf HG01565 HG01566 HG01571 HG01572 HG01577 HG01578 HG01892 HG01893 HG01917 HG01918 HG01920 HG01921 HG01923 HG01924 HG01926 HG01927 HG01932 HG01933 HG01935 HG01938 HG01939 HG01941 HG01942 HG01944 HG01945 HG01947 HG01948 HG01950 HG01951 HG01953 HG01954 HG01961 HG01965 HG01967 HG01968 HG01970 HG01971 HG01973 HG01974 HG01976 HG01977 HG01979 HG01980 HG01982 HG01991 HG01992 HG01995 HG01997 HG02002 HG02003 HG02008 HG02089 HG02090 HG02104 HG02105 HG02146 HG02147 HG02252 HG02253 HG02259 > Data/PEL.$GENE.vcf
sed -i -e 's/\t\./\t0|0/g' Data/PEL.$GENE.vcf

$VCFLIB/vcfkeepsamples Data/ALL.$GENE.vcf NA18524 NA18525 NA18526 NA18527 NA18528 NA18529 NA18530 NA18531 NA18532 NA18533 NA18534 NA18535 NA18536 NA18537 NA18538 NA18539 NA18540 NA18541 NA18542 NA18543 NA18544 NA18545 NA18546 NA18547 NA18548 NA18549 NA18550 NA18552 NA18553 NA18555 NA18557 NA18558 NA18559 NA18560 NA18561 NA18562 NA18563 NA18564 NA18565 NA18566 NA18567 NA18569 NA18570 NA18571 NA18572 NA18573 NA18574 NA18575 NA18576 NA18577 NA18579 NA18580 NA18582 NA18583 NA18591 NA18592 NA18593 NA18595 NA18596 NA18597 > Data/CHB.$GENE.vcf
sed -i -e 's/\t\./\t0|0/g' Data/CHB.$GENE.vcf

$VCFLIB/vcfkeepsamples Data/ALL.$GENE.vcf HG01112 HG01113 HG01119 HG01121 HG01122 HG01124 HG01125 HG01130 HG01131 HG01133 HG01134 HG01136 HG01137 HG01139 HG01140 HG01142 HG01148 HG01149 HG01250 HG01251 HG01253 HG01254 HG01256 HG01257 HG01259 HG01260 HG01269 HG01271 HG01272 HG01274 HG01275 HG01277 HG01278 HG01280 HG01281 HG01284 HG01341 HG01342 HG01344 HG01345 HG01347 HG01348 HG01350 HG01351 HG01353 HG01354 HG01356 HG01357 HG01359 HG01360 HG01362 HG01363 HG01365 HG01366 HG01369 HG01372 HG01374 HG01375 HG01377 HG01378 > Data/CLM.$GENE.vcf
sed -i -e 's/\t\./\t0|0/g' Data/CLM.$GENE.vcf

# clean up
rm $VCF $VCF.tbi ALL.chr$CHROM.tmp.vcf ALL.$GENE.tmp.vcf


echo Downloading recombination map
mkdir Map
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/PEL_omni_recombination_20130507.tar &> /dev/null
tar -C Map -xvf PEL_omni_recombination_20130507.tar &> /dev/null
zcat Map/PEL/PEL-$CHROM-final.txt.gz > Data/genetic_map_chrom$CHROM.map
rm -rf Map
rm PEL_omni_recombination_20130507.tar


echo Done!
ls -lh Data/* > Data/download.log
echo Open Data/download.log to see what files have been generated.





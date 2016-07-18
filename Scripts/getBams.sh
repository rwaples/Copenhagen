
echo Retrieving file names...

# EDAR is our target gene
# chr2:109,510,927-109,605,828

# sample size
NS=20

# chrom start end
CHROM=2
START=108000000
END=1110000000

# get IDs 
for POP in LWK TSI CHB PEL;
do

	INPUT=$POP.txt
	INDS=`cat $INPUT`
	OUTPUT=$POP.BAMs.txt
	> $OUTPUT
	> tmp

	for i in $INDS;
	do
		grep $i Files/phase3_bamlist.txt >> tmp
	done
	# wc -l tmp
	head -n $NS tmp > $OUTPUT 
	rm tmp
done

echo Downloading BAM files...
# download and index bams
for POP in LWK TSI CHB PEL;
do
	mkdir Data/$POP.BAMs
	echo $POP
	INDLIST=`cat $POP.BAMs.txt`
	for i in $INDLIST;
	do
		NAME=`echo -n $i | tail -c 58`
		echo $NAME
		samtools view -s 0.25 -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/$i $CHROM:$START-$END > Data/$POP.BAMs$NAME 2> /dev/null
		#samtools index $POP.BAMs/$NAME
	done
done

echo Removing index files...
rm *.bai

exit




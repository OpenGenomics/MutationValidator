#!/bin/sh

. /broad/tools/scripts/useuse



use .zlib-1.2.6 #having to make this wrapper and hack a new python version because firehose needed this specific use statement and then optional arguments became a pain. :(
use Samtools
#use .samtools-0.1.9

#unuse .python-2.7.1-sqlite3-rtrees
#reuse .python-3.4.2

LIBDIR=$1
MAFSNP=$2
MAFINDEL=$3
WEXTUMOR=$4
WEXNORMAL=$5
WGSTUMOR=$6
WGSNORMAL=$7
RNATUMOR=$8
TARGT=$9
TARGN=$10
LPT=$11 #lowpass
LPN=$12
OT=$13 #other
ON=$14
PAIRID=$15

RNATYPE="None"

#check the alignment of the rna
if [ -f $RNATUMOR ]; then
	echo "rna tumor file exists"
	CHECK=$(samtools view -H $RNATUMOR | grep SN:chr1)

	if [ -z $CHECK ]; then 
		echo "file is aligned correctly"
		RNATYPE="hg19"
	else 
		echo "file is aligned with chr in front"
		RNATYPE="hg19-chr"
	fi
fi

echo "python $LIBDIR/validation_wrapper_firehose_library_hack.py --mafsnp $MAFSNP --mafindel $MAFINDEL --wextumor $WEXTUMOR --wexnormal $WEXNORMAL --wgstumor $WGSTUMOR --wgsnormal $WGSNORMAL --rnatumor $RNATUMOR --targetedtumor $TARGT --targetednormal $TARGN --lowpasstumor $LPT --lowpassnormal $LPN --othertumor $OT --othernormal $ON --out $PAIRID --rnatype $RNATYPE"

python $LIBDIR/validation_wrapper_firehose_library_hack.py --mafsnp $MAFSNP --mafindel $MAFINDEL --wextumor $WEXTUMOR --wexnormal $WEXNORMAL --wgstumor $WGSTUMOR --wgsnormal $WGSNORMAL --rnatumor $RNATUMOR --targetedtumor $TARGT --targetednormal $TARGN --lowpasstumor $LPT --lowpassnormal $LPN --othertumor $OT --othernormal $ON --out $PAIRID --rnatype $RNATYPE


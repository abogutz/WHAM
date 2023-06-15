#!/bin/bash




#  _       ____  _____    __  _____
# | |     / / / / /   |  /  |/  / /
# | | /| / / /_/ / /| | / /|_/ / / 
# | |/ |/ / __  / ___ |/ /  / /_/  
# |__/|__/_/ /_/_/  |_/_/  /_(_)
#
# WHAM: WGBS Heterogeneity Analysis at the Molecular level
# 3-part toolkit for analysis of DNAme distribution
# 1: Lollipop visualization of by-read methylation calls
# 2: Heatmap of read-level methylation distribution
# 3: Modality test for non-unimodal distribution
# Requires samtools, bedtools, R, bedToBigBed, bedGraphToBigWig
# Within R: requires diptest, foreach, and doParallel




SCRIPTS_DIR="project/def-mlorincz/scripts/misc/WGBS/"
LOLLY_SCRIPT=$SCRIPTS_DIR"lolly.awk"
DIP_AWK_SCRIPT=$SCRIPTS_DIR"diptest-bin.awk"
HEAT_AWK_SCRIPT=$SCRIPTS_DIR"heatmap-bin.awk"
BIGLOLLY_AS=$SCRIPTS_DIR"bigLolly-size.as"
R_SCRIPT=$SCRIPTS_DIR"DipTest-parr.R"

CHR_SIZES="/mnt/d/Data/Annotations/mm10/mm10.chrom.sizes"


SCRATCH_DIR=""
TEMP1=$SCRATCH_DIR"/temp"
TEMP2=$SCRATCH_DIR"/temp2"

THREADS=$SLURM_CPUS_PER_TASK

MAPQ=40
MIN_CPG=4
DIPTEST_BINSIZE=100
HEAT_GENOME_BINSIZE=25
METH_BINS=5
MAX_READS=10
COLOR_BINS=5

INPUT=$1
NAME=$(basename $INPUT .bam)
LOLLY_OUTPUT=$NAME"_lolly.bb"
DIPTEST_OUTPUT=$NAME"-"$GENOME_BINSIZE"bp"$MIN_CPG"CpG-diptest.bw"

GENOME_BED=$SCRATCH_DIR"/"$(basename $CHR_SIZES .sizes)".bounds"
awk 'OFS="\t"{print $1, 0, $2}' $CHR_SIZES > $TEMP1
sort -k1,1 -k2,2n $TEMP1 > $GENOME_BED

function initializeHub () {
	HUB="Track_Hub/"
	GENOME_DIR=$HUB"mm10/"
	TRACKDB=$GENOME_DIR"trackDb.txt"
	mkdir $HUB
	mkdir $GENOME_DIR
	printf "hub <HubNameWithoutSpace>\nshortLabel <max 17 char, display on side>\nlongLabel Hub to display <fill> data at UCSC\ngenomesFile genomes.txt\nemail <email-optional>" > $HUB/hub.txt
	printf "genome mm10\ntrackDb mm10/trackDb.txt" > $HUB/genomes.txt
}

function makeLollies () {
	echo "Starting Lolly Generation"
	samtools view -q $MAPQ $INPUT | awk -f $LOLLY_SCRIPT | sort -k1,1 -k2,2n -T $SCRATCH_DIR > $TEMP1
	bedToBigBed -as=$BIGLOLLY_AS -type=bed9+1 $TEMP1 $CHR_SIZES $LOLLY_OUTPUT
	rm $TEMP1
	printf "track %s\nshortLabel %s\nlongLabel %s\ntype bigLolly\nbigDataUrl %s\nvisibility full\nlollyNoStems on\nlollySizeField 10\nmaxWindowToDraw 20000\n\n" $LOLLY_OUTPUT $LOLLY_OUTPUT $LOLLY_OUTPUT $LOLLY_OUTPUT | tee -a $TRACKDB
	echo "Finished Lolly Generation"
}


function dipTest () {
	echo "Starting Diptest Calculations"
	REF_BED=$SCRATCH_DIR"/ref"-$DIPTEST_BINSIZE"bp.bed"
	
	echo "Making genomic windows..."
	bedtools makewindows -w $DIPTEST_BINSIZE -b $GENOME_BED > $REF_BED

	echo "Counting Bins..."
	samtools view -q $MAPQ $INPUT | awk -f $DIP_AWK_SCRIPT -v thresh=$MIN_CPG > $TEMP1

	echo "Mapping..."
	bedtools map -a $REF_BED -b $TEMP1 -c 4 -o collapse | awk 'OFS="\t"{if($4 != ".") {print $0}}' > $TEMP2

	# R CODE STUFF
	echo "Calculating Dip p-values..."
	Rscript $R_SCRIPT $TEMP2 $TEMP1 $THREADS
	sed -i 's/\"//g' $TEMP1
	echo "Generating bigwig..."
	bedGraphToBigWig $TEMP1 $CHR_SIZES $DIPTEST_OUTPUT
	rm $REF_BED $TEMP1 $TEMP2
	printf "track %s\nshortLabel %s\nlongLabel %s\ntype bigWig\nbigDataUrl %s\ncolor 0,0,0\nvisibility full\nmaxHeightPixels 100:60:25\nautoScale on\nalwaysZero on\nyLineOnOff on\nyLineMark 1.3\n\n" $DIPTEST_OUTPUT $DIPTEST_OUTPUT $DIPTEST_OUTPUT $DIPTEST_OUTPUT | tee -a $TRACKDB
}


function heatmap () {
	let METH_BINS_SIZE=100/$METH_BINS
	let COLOR_BINS_SIZE=$MAX_READS/$COLOR_BINS
	PRIORITY=1
	REF_BED="ref"-$HEAT_GENOME_BINSIZE"bp.bed"
	
	printf "track %s\ncontainer multiWig\nshortLabel %s\nlongLabel %s\ntype bigWig\nvisibility full\nmaxHeightPixels 100:60:25\nconfigurable on\nviewLimits 0:100\nalwaysZero on\naggregate solidOverlay\nshowSubtrackColorOnUi on\npriority 1.0\n\n" $NAME $NAME $NAME | tee -a $TRACKDB

	echo "Making genomic windows..."
	bedtools makewindows -w $HEAT_GENOME_BINSIZE -b $GENOME_BED > $REF_BED

	echo "Counting Bins..." # First level of binning - put reads into methylation bins
	samtools view -q $MAPQ $INPUT | awk -f $HEAT_AWK_SCRIPT -v bins=$METH_BINS -v thresh=$MIN_CPG

	for (( BIN=100; BIN>=0; BIN-=$METH_BINS_SIZE ))
	do
		echo "Processing bin "$BIN
		FILE="bin"$BIN
		BEDGRAPH=$SCRATCH_DIR"/"$FILE".bedgraph"
		
		echo "Mapping..." # Second binning - pileup reads in each genomic window
		bedtools map -a $REF_BED -b $FILE -c 1 -o count > $BEDGRAPH
		rm $FILE
		
		echo "Binning by depth..." # Third binning - determine color bin from read depth
		awk -v bin=$BIN -v colorBinSize=$COLOR_BINS_SIZE -v maxReads=$MAX_READS 'OFS="\t"{
			done=0;
			if($4 > 0) {
				for ( x=colorBinSize; x<maxReads; x+=colorBinSize ) {
					if ($4 < x) {
						print $1, $2, $3, bin > "bin"bin"-"x".bedgraph";
						done=1;
						break;
					}
				}
				if (!done) {
					print $1, $2, $3, bin > "bin"bin"-"maxReads".bedgraph";
				}
			}
		}' $BEDGRAPH
		rm $BEDGRAPH
		
		# Generate entire genome bins for 0-0
		awk -v bin=$BIN 'OFS="\t"{
			print $0, bin;
		}' $GENOME_BED > $FILE"-0.bedgraph"
		
		# Generate trackdb track layout and coloring 
		for BG in $SCRATCH_DIR"/"$FILE-*.bedgraph #(( CURR_BIN=0; CURR_BIN<=$MAX_READS; CURR_BIN+=$COLOR_BINS_SIZE ))
		do
			CURR_BIN=${BG//$FILE-/}
			CURR_BIN=${CURR_BIN//.bedgraph/}
			BW_NAME=${BG//.bedgraph/.bw}
			BW=$GENOME_DIR$BW_NAME
			if [[ -f $BG ]]; then
				bedGraphToBigWig $BG $CHR_SIZES $BW
				rm $BG
				R=$(echo "scale=5;255*(3-(3*$CURR_BIN/$MAX_READS))" | bc)
				R=$(echo "scale=0;$R/1" | bc)
				R=$(( 255 < $R ? 255 : $R ))
				G=$(echo "scale=5;255*(2-(3*$CURR_BIN/$MAX_READS))" | bc)
				G=$(echo "scale=0;$G/1" | bc)
				G=$(( 0 > $G ? 0 : $G ))
				G=$(( 255 < $G ? 255 : $G ))
				B=$(echo "scale=5;255*(1-(3*$CURR_BIN/$MAX_READS))" | bc)
				B=$(echo "scale=0;$B/1" | bc)
				B=$(( 0 > $B ? 0 : $B ))
				COLOR=$R","$G","$B
				printf "\ttrack %s\n\tparent %s\n\tshortLabel %s\n\tlongLabel %s\n\ttype bigWig\n\tbigDataUrl %s\n\tcolor %s\n\tpriority %s\n\n" $BW_NAME $NAME $BW_NAME $BW_NAME $BW_NAME $COLOR $PRIORITY | tee -a $TRACKDB
				((PRIORITY++))
			fi
		done
	done
}






#rm -r $SCRATCH_DIR


#! /bin/bash
#SBATCH --account=def-mlorincz           # required (format def-name)
#SBATCH --cpus-per-task=10                        # number of cpus
#SBATCH --mem-per-cpu=4G                 # memory; default unit is megabytes
#SBATCH --time=01-12:00                   # time (DD-HH:MM)
#SBATCH --mail-user=aaron.bogutz@ubc.ca
#SBATCH --mail-type=ALL

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
# If creating traditional plots, additionally requires bismark and Deeptools bamCoverage

## Scripts Locations - CHANGE TO ACTUAL LOCATION ##
SCRIPTS_DIR="/project/def-mlorincz/scripts/misc/WHAM/"

PE_PARSER=$SCRIPTS_DIR"PE-methParser.awk"
LOLLY_SCRIPT=$SCRIPTS_DIR"lolly.awk"
DIP_AWK_SCRIPT=$SCRIPTS_DIR"diptest-bin.awk"
HEAT_AWK_SCRIPT=$SCRIPTS_DIR"heatmap-bin.awk"
BIGLOLLY_AS=$SCRIPTS_DIR"bigLolly-size.as"
R_SCRIPT=$SCRIPTS_DIR"DipTest.R"
CONFIG=$SCRIPTS_DIR"ComputeCanada.config"
source $CONFIG

BAMCOVERAGE="/project/def-mlorincz/scripts/utilities/miniconda3/bin/bamCoverage"



## Default Values ##
MAPQ=40
MIN_CPG=4
DIPTEST_BINSIZE=100
HEAT_GENOME_BINSIZE=25
METH_BINS=5
MAX_READS=10
COLOR_BINS=5
SCRATCH_DIR="./"
TEMP1=$SCRATCH_DIR"/temp"
TEMP2=$SCRATCH_DIR"/temp2"
PE_BAM=$SCRATCH_DIR"/temp.bam"
THREADS=$SLURM_CPUS_PER_TASK
CHR_SIZES="/project/def-mlorincz/reference_genomes/mm10/mm10.sizes"
TRAD=0
LOLLY=1
HEATMAP=1

## Help Messages ##
HELP="USAGE:\t $(basename $0) [OPTIONS] -h for help"
HELP_FULL="\n$HELP\n
\nThis set of scripts will perform 3 analyses of Bismark-aligned data:\n\t1: Lollipop visualization of by-read methylation calls\n\t2: Heatmap of read-level methylation distribution\n\t3: Modality test for non-unimodal distribution\nRequires samtools, bedtools, bedToBigBed, bedGraphToBigWig, R.\nWithin R: requires diptest, foreach, and doParallel\n\n
OPTIONS:\n\t
-h\tPrints help page.\n\t
-d\tCheck dependencies and exit.\n\t
-i\tInput file. Must be a Bismark aligned .bam file. REQUIRED\n\t
-s\tScratch directory. Default=./\n\t
-t\tNumber of threads to use. Default=SLURM_THREADS\n\t
-q\tMinimum mapping quality for reads. Default=$MAPQ\n\t
-C\tMinimum CpGs for reads. Default=$MIN_CPG\n\t
-D\tGenome bin size for Diptest. Default=$DIPTEST_BINSIZE\n\t
-H\tGenome bin size for Heatmp. Default=$HEAT_GENOME_BINSIZE\n\t
-B\tNumber of methylation bins for Heatmap. Default=$METH_BINS\n\t
-R\tTop end of teads for color in Heatmap. Default=$MAX_READS\n\t
-c\tNumber of color bins for Heatmap. Default=$COLOR_BINS\n\t
-z\tChromosome sizes file. Default= $CHR_SIZES\n\t
-p\tCreate traditional methylation and coverage tracks. Default=OFF\n\t
-l\tDon't create lollipop tracks. Default=ON\n\t
-m\tDon't create heatmap track. Default=ON\n\t
-b\tReference Bed File. Will ONLY perform diptest."


OPTIONS="hi:q:C:D:H:B:R:c:s:t:z:dplmb:"

function parseOptions () {
	if ( ! getopts $OPTIONS opt); then
		echo -e $HELP
		exit 1
	fi

	while getopts $OPTIONS opt; do
		case $opt in
			h) #open the HELP menu
				echo -e $HELP_FULL | fold -s
				exit
				;;
			i) #set input file
				INPUT=${OPTARG}
				NAME=$(basename $INPUT .bam)
				;;
			q) #minimum MapQ
				MAPQ=${OPTARG}
				;;
			C) #minimum CpGs in read
				MIN_CPG=${OPTARG}
				;;
			D) #binsize for Diptest
				DIPTEST_BINSIZE=${OPTARG}
				;;
			H) #binsize for heatmap
				HEAT_GENOME_BINSIZE=${OPTARG}
				;;
			B) #number of methylation bins for heatmap
				METH_BINS=${OPTARG}
				;;
			R) #max read depth for heatmap
				MAX_READS=${OPTARG}
				;;
			c) #number of colour bins
				COLOR_BINS=${OPTARG}
				;;
			s) #scratch directory
				SCRATCH_DIR=${OPTARG}
				mkdir $SCRATCH_DIR
				TEMP1=$SCRATCH_DIR"/temp"
				TEMP2=$SCRATCH_DIR"/temp2"
				PE_BAM=$SCRATCH_DIR"/temp.bam"
				;;
			t) #threads
				THREADS=${OPTARG}
				;;
			z) #chromosome sizes file
				CHR_SIZES=${OPTARG}
				;;
			d) #check Dependencies
				checkDependencies
				exit 0;
				;;
			p) #make traditional plots
				TRAD=1
				;;
			l) #Don't make lollipops
				LOLLY=0
				;;
			m) #don't make heatmap
				HEATMAP=0
				;;
			b) #input bed file
				REF_BED=${OPTARG}
				HEATMAP=0
				LOLLY=0
				;;
			\?)
				echo -e "\n###############\nERROR: Invalid Option! \nTry '$(basename $0) -h' for help.\n###############" >&2
				exit 1
				;;
		esac
	done
	if [[ ! -f $INPUT ]]; then
		echo "Not a valid input file."
		exit 1
	fi
}

function checkDependencies () {
	DEPENDENCIES=(bedtools awk samtools Rscript bedToBigBed bedGraphToBigWig)
	if [[ $TRAD == 1 ]] ; then #Create Traditional Plots
		DEPENDENCIES=(${DEPENDENCIES[@]} bismark $BAMCOVERAGE)
	fi
	echo -e "Checking Dependencies:"
	EXIT=0
	for COMMAND in "${DEPENDENCIES[@]}"; do
		echo -e "[checkDependencies] $COMMAND..."
		command -v $COMMAND > /dev/null 2>&1 || {
			echo -e >&2 "\t\t$COMMAND not found!"
			EXIT=1
		}
	done
	if [[ $EXIT = 1 ]] ; then
		exit 1
	fi
}

function initializeHub () {
	HUB="Track_Hub/"
	GENOME=$(basename $CHR_SIZES | cut -d. -f1)
	GENOME_DIR=$HUB$GENOME"/"
	TRACKDB=$GENOME_DIR"trackDb.txt"
	mkdir $HUB
	mkdir $GENOME_DIR
	printf "hub <HubNameWithoutSpace>\nshortLabel <max 17 char, display on side>\nlongLabel Hub to display <fill> data at UCSC\ngenomesFile genomes.txt\nemail <email-optional>" > $HUB/hub.txt
	printf "genome %s\ntrackDb %s/trackDb.txt" $GENOME $GENOME > $HUB/genomes.txt
	
	LOLLY_OUTPUT=$GENOME_DIR$NAME"_lolly.bb"
	DIPTEST_OUTPUT=$GENOME_DIR$NAME"-"$DIPTEST_BINSIZE"bp"$MIN_CPG"CpG-diptest.bw"
	METHYL_OUTPUT=$GENOME_DIR$NAME"_methylation.bw"
	COVERAGE_OUTPUT=$GENOME_DIR$NAME"_coverage.bw"
}

function makeLollies () {
	echo "Starting Lolly Generation"
	echo "Converting to lolly..."
	samtools view -q $MAPQ $INPUT | awk -f $LOLLY_SCRIPT > $TEMP2
	sort -k1,1 -k2,2n --parallel=$THREADS --buffer-size="3G" -T $SCRATCH_DIR $TEMP2 > $TEMP1
	echo "Compressing lollies..."
	bedToBigBed -as=$BIGLOLLY_AS -type=bed9+1 $TEMP1 $CHR_SIZES $LOLLY_OUTPUT
#	rm $TEMP1
	LOLLY_NAME=$(basename $LOLLY_OUTPUT)
	printf "track %s\nshortLabel %s\nlongLabel %s\ntype bigLolly\nbigDataUrl %s\nvisibility full\nautoScale on\nlollyNoStems on\nlollySizeField 10\nmaxWindowToDraw 20000\n\n" $LOLLY_NAME $LOLLY_NAME $LOLLY_NAME $LOLLY_NAME | tee -a $TRACKDB
	echo "Finished Lolly Generation"
}

function dipTest () {
	echo "Starting Diptest Calculations"
	
	echo "Making genomic windows..."
	if [[ $REF_BED == "" ]] ; then # If not using input bed, then bin genome
		REF_BED=$SCRATCH_DIR"/ref"-$DIPTEST_BINSIZE"bp.bed"
		bedtools makewindows -w $DIPTEST_BINSIZE -b $GENOME_BED > $REF_BED
	else # Using input bed file
		REF_BED_NAME=${REF_BED%%.*}
		DIPTEST_OUTPUT=$GENOME_DIR$NAME"_"$REF_BED_NAME"-"$MIN_CPG"CpG-diptest.bw"
	fi

	echo "Counting Bins..."
	samtools view -q $MAPQ $INPUT | awk -f $DIP_AWK_SCRIPT -v thresh=$MIN_CPG > $TEMP1

	echo "Mapping..."
	bedtools map -a $REF_BED -b $TEMP1 -c 4 -o collapse | awk 'OFS="\t"{if($4 != ".") {print $0}}' > $TEMP2

	# R CODE STUFF
	echo "Calculating Dip p-values..."
	Rscript $R_SCRIPT $TEMP2 $TEMP1 $THREADS
	sed -i 's/\"//g' $TEMP1
	echo "Generating bigwig..."
	sort -k1,1 -k2,2n $TEMP1 > $TEMP2
	bedGraphToBigWig $TEMP2 $CHR_SIZES $DIPTEST_OUTPUT
	rm $TEMP1 $TEMP2
	DIPTEST_NAME=$(basename $DIPTEST_OUTPUT)
	printf "track %s\nshortLabel %s\nlongLabel %s\ntype bigWig\nbigDataUrl %s\ncolor 255,0,0\nvisibility full\nmaxHeightPixels 100:60:25\nautoScale on\nalwaysZero on\nyLineOnOff on\nyLineMark 1.3\n\n" $DIPTEST_NAME $DIPTEST_NAME $DIPTEST_NAME $DIPTEST_NAME | tee -a $TRACKDB
}

function heatmap () {
	let METH_BINS_SIZE=100/$METH_BINS
	let COLOR_BINS_SIZE=$MAX_READS/$COLOR_BINS
	PRIORITY=1
	REF_BED=$SCRATCH_DIR"ref"-$HEAT_GENOME_BINSIZE"bp.bed"
	BW_DIR=$GENOME_DIR$NAME"/"
	mkdir $BW_DIR
	
	printf "track %s\ncontainer multiWig\nshortLabel %s\nlongLabel %s\ntype bigWig\nvisibility full\nmaxHeightPixels 100:60:25\nconfigurable on\nviewLimits 0:100\nalwaysZero on\naggregate solidOverlay\nshowSubtrackColorOnUi on\npriority 1.0\n\n" $NAME $NAME $NAME | tee -a $TRACKDB

	echo "Making genomic windows..."
	bedtools makewindows -w $HEAT_GENOME_BINSIZE -b $GENOME_BED > $REF_BED

	echo "Counting Bins..." # First level of binning - put reads into methylation bins
	samtools view -q $MAPQ $INPUT | awk -f $HEAT_AWK_SCRIPT -v bins=$METH_BINS -v thresh=$MIN_CPG # TODO should this be done in scratch?

	for (( BIN=100; BIN>=0; BIN-=$METH_BINS_SIZE ))
	do
		echo "Processing bin "$BIN
		FILE="bin"$BIN
		BEDGRAPH=$SCRATCH_DIR"/"$FILE".bedgraph"
		
		echo "Mapping..." # Second binning - pileup reads in each genomic window
		bedtools map -a $REF_BED -b $FILE -c 1 -o count > $BEDGRAPH
		rm $FILE
		
		echo "Binning by depth..." # Third binning - determine color bin from read depth
		pushd $SCRATCH_DIR # Move to scratch directory for bin bigwig creations
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
		popd # Return from scratch directory
		
		# Generate trackdb track layout and coloring 
		for BG in $SCRATCH_DIR"/"$FILE-*.bedgraph #(( CURR_BIN=0; CURR_BIN<=$MAX_READS; CURR_BIN+=$COLOR_BINS_SIZE ))
		do
			CURR_BIN=$(basename $BG .bedgraph)
			CURR_BIN=${CURR_BIN//$FILE-/}
			BW_NAME=$(basename $BG .bedgraph).bw
			BW=$BW_DIR$BW_NAME
			if [[ -f $BG ]]; then
				sort -k1,1 -k2,2n $BG > $TEMP1
				bedGraphToBigWig $TEMP1 $CHR_SIZES $BW
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
				printf "\ttrack %s\n\tparent %s\n\tshortLabel %s\n\tlongLabel %s\n\ttype bigWig\n\tbigDataUrl %s\n\tcolor %s\n\tpriority %s\n\n" $BW_NAME $NAME $BW_NAME $BW_NAME $NAME"/"$BW_NAME $COLOR $PRIORITY | tee -a $TRACKDB
				((PRIORITY++))
			fi
		done
	done
}

function makeBounds () { # Take chr sizes file and make bed file
	GENOME_BED=$SCRATCH_DIR"/"$(basename $CHR_SIZES .sizes)".bounds"
	awk 'OFS="\t"{print $1, 0, $2}' $CHR_SIZES > $TEMP1
	sort -k1,1 -k2,2n $TEMP1 > $GENOME_BED
}

function parsePE () { # Combine Methylation strings from PE reads into a single entry
	FLAG=$(samtools view $INPUT | head -n 1 | cut -f2)
	if [[ $((FLAG&1)) == 1 ]] ; then #Paired-end
		PAIRED=true
	else
		PAIRED=false
	fi
	echo "Data are Paired-end:" $PAIRED
	if [[ $PAIRED == true ]] ; then
		if [[ $REF_BED != "" ]] ; then # If using a reference bed file, only keep reads within 1kb of regions
			bedtools slop -i $REF_BED -g $CHR_SIZES -b 1000 > $TEMP1
			samtools view -bh -q $MAPQ -L $TEMP1 $INPUT > $PE_BAM
			INPUT=$PE_BAM
		fi
		echo "Sorting PE BAM by read name..."
		samtools sort -@ $THREADS -m 3G -T $SCRATCH_DIR -n -o $TEMP1 $INPUT
		echo "Combining Methylation for PE reads..."
		samtools view -q $MAPQ $TEMP1 | awk -f $PE_PARSER > $TEMP2
		samtools view -H $INPUT > $TEMP1
		echo "Compressing PE Data..."
		cat $TEMP1 $TEMP2 | samtools view -@ $THREADS -bh -o $TEMP1
		echo "Sorting PE Data..."
		samtools sort -@ $THREADS -m 3G -T $SCRATCH_DIR -o $PE_BAM $TEMP1
		INPUT=$PE_BAM
	fi
}

function makeTradPlots () {
	samtools index $INPUT
	$BAMCOVERAGE --minMappingQuality $MAPQ --outFileFormat bigwig -p $THREADS -b $INPUT -o $COVERAGE_OUTPUT
	bismark_methylation_extractor -o $SCRATCH_DIR --gzip --multicore $THREADS --bedGraph --mbias_off $INPUT
	zcat $SCRATCH_DIR$NAME".bedGraph.gz" > $TEMP1
	tail -n +2 $TEMP1 > $TEMP2
	sort -k1,1 -k2,2n $TEMP2 > $TEMP1
	bedGraphToBigWig $TEMP1 $CHR_SIZES $METHYL_OUTPUT
	METHYL_NAME=$(basename $METHYL_OUTPUT)
	printf "track %s\nshortLabel %s\nlongLabel %s\ntype bigWig\nbigDataUrl %s\ncolor 0,0,0\nvisibility full\nmaxHeightPixels 100:60:25\nautoScale on\nalwaysZero on\nyLineOnOff on\nyLineMark 1.3\n\n" $METHYL_NAME $METHYL_NAME $METHYL_NAME $METHYL_NAME | tee -a $TRACKDB
	COVERAGE_NAME=$(basename $COVERAGE_OUTPUT)
	printf "track %s\nshortLabel %s\nlongLabel %s\ntype bigWig\nbigDataUrl %s\ncolor 0,0,0\nvisibility full\nmaxHeightPixels 100:60:25\nautoScale on\nalwaysZero on\nyLineOnOff on\nyLineMark 1.3\n\n" $COVERAGE_NAME $COVERAGE_NAME $COVERAGE_NAME $COVERAGE_NAME | tee -a $TRACKDB
}

## Actually Run Stuff ##
loadModules
parseOptions $@
makeBounds
initializeHub
if [[ $TRAD == 1 ]] ; then #Create Traditional Plots
	makeTradPlots
fi
parsePE
if [[ $LOLLY == 1 ]] ; then #Create Lollipops
	makeLollies
fi
if [[ $HEATMAP == 1 ]] ; then #Create Heatmap
	heatmap
fi
dipTest

rm -r $SCRATCH_DIR


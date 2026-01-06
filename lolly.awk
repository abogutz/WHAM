#! /bin/awk

# Helper script to convert Bismark BAM files into bigLolly files for UCSC
# Made by Aaron
# Last Updated 2023-05-10

BEGIN {
	OFS="\t";
	black="0,0,0";
	white="255,255,255";
	green="0,255,0";
	red="255,0,0"
	grey="150,150,150";
	XM=0;
	if( ! size ) { size=2; } # Default size value - less than 2 is somewhat unreadable
	if( ! capSize ) { capSize=2; } # Size for start and end flags
	if( ! minDist ) { minDist=20; } # Arbitrary distance between reads on same "strand"
	if( ! gapStep ) { gapStep=5; } # Arbitrary distance gap indicators
	if( ! gapSize ) { gapSize=1; } # Make gaps itty bitty
}{
	while ( $3 !~ /chr[0-9XY]*$/ && done>0 ) { # Ignore non-canonical chr
		done=getline;
	}
	if(XM == 0) { # Locate XM string
		for(x=12; x<20; x++) { #Locate XM Tag
			if($x ~ /^XM:Z:/) {
				XM=x;
				break;
			}
		}
		if (XM == 0) {
			print "Failed to locate Tag" > "/dev/stderr";
			exit 1;
		}
	}
	if ( $3 != chr ) { # New chromosome - reset strand count
		chr=$3;
		for(x=1; x<1000; x++) {
			end[x]=0;
		}
	}
	methCalls=substr($XM, 6) # Remove leading characters in meth calls - which should always be in column 16
	calls=split(methCalls, len, /[zZ]/, meth); # Split methylation calls into inter-CpG distances (len) and CpG meth calls (meth)
	if ( calls > 1 ) {
		start=$4;
		for(strand=1; strand<1000; strand++) {
			if( start>end[strand]+minDist ) {
				break;
			}
		}
		if(strand < 1000 ) {
			print chr, start-1, start, "start", strand, ".", 0, 0, green, capSize; # print start in green
			for ( x=1 ; x<calls ; x++ ) {
				start+=length(len[x]);
				if ( meth[x] == "z" ) { # unmethylated
					colour = white;
				} else { # methylated
					colour = black;
				}
				print chr, start-1, start, meth[x], strand, ".", 0, 0, colour, size;
				start++; # Increment over Cytosine
			}
			start+=length(len[x]); # Add final interval
			end[strand]=start;
			print chr, start-1, start, "end", strand, ".", 0, 0, red, capSize; # print end in red
			gap=match(methCalls, /\++/);
			if(gap != 0) {
				for(start=0; start<RLENGTH; start+=gapStep) {
								print chr, $4+gap+start-1, $4+gap+start, "gap", strand, ".", 0, 0, grey, gapSize; # gap in grey
				}
			}
		}
	}
}
#! /bin/awk

# Helper script to convert Bismark BAM files into bigLolly files for UCSC
# Made by Aaron
# Last Updated 2026-07-07 using Marlet's lolly variant (light-blue caps and per-bp markers)

BEGIN {
	OFS="\t";
	black="0,0,0";
	white="255,255,255";
	start_col="210,230,255"; # light-blue start cap
	end_col="210,230,255"; # light-blue end cap
	grey="150,150,150";
	XM=0;
	XG=0;
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
	if(XG == 0) { 
		for(x=12; x<20; x++) {
			if($x ~ /^XG:Z:/) {
				XG=x;
				break;
			}
		}
		if (XG == 0) {
			print "Failed to locate XG Tag" > "/dev/stderr";
			exit 1;
		}
	}
	if ( $3 != chr ) { # New chromosome - reset strand count
		chr=$3;
		for(x=1; x<1000; x++) {
			end[x]=0;
		}
	}
	methCalls=substr($XM, 6) # Remove leading characters in meth calls - which should always be in column 16 (BUT IT ISNT)
	calls=split(methCalls, len, /[zZ]/, meth); # Split methylation calls into inter-CpG distances (len) and CpG meth calls (meth)
	if ( calls > 1 ) {
		start=$4;
		for(strand=1; strand<1000; strand++) {
			if( start>end[strand]+minDist ) {
				break;
			}
		}
		if(strand < 1000 ) {
			# Offset start by the genome strand the read aligned to (Marlet)
			if($XG == "XG:Z:CT") {
				start=$4;
			} else { # XG:Z:GA
				start=$4-1;
			}
			print chr, start-2, start-1, "start", strand, ".", 0, 0, start_col, capSize; # light-blue start cap
			for ( x=1 ; x<calls ; x++ ) {
				num_bp=length(len[x]); # per-bp markers across the inter-CpG interval (Marlet)
				for ( i=1 ; i<num_bp ; i+=3 ) {
					print chr, start+i-1, start+i, "bp", strand, ".", 0, 0, grey, 0;
				}
				start+=length(len[x]);
				if ( meth[x] == "z" ) { # unmethylated
					colour = white;
				} else { # methylated
					colour = black;
				}
				print chr, start-1, start, meth[x], strand, ".", 0, 0, colour, size;
				start++; # Increment over Cytosine
			}
			num_bp=length(len[x]); # per-bp markers for the final interval (Marlet)
			for ( i=1 ; i<num_bp ; i+=3 ) {
				print chr, start+i-1, start+i, "bp", strand, ".", 0, 0, grey, 0;
			}
			start+=length(len[x]); # Add final interval
			end[strand]=start;
			print chr, start, start+1, "end", strand, ".", 0, 0, end_col, capSize; # light-blue end cap
			gap=match(methCalls, /\++/);
			if(gap != 0) {
				for(start=0; start<RLENGTH; start+=gapStep) {
								print chr, $4+gap+start-1, $4+gap+start, "gap", strand, ".", 0, 0, grey, gapSize; # gap in grey
				}
			}
		}
	}
}

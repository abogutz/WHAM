#! /bin/awk

# Helper script to convert Bismark BAM files into bigLolly files for UCSC
# Made by Aaron
# Last Updated 2023-05-10

BEGIN{
	OFS="\t";
	chr="";
	black="0,0,0";
	white="255,255,255";
	green="0,255,0";
	red="255,0,0"
	grey="100,100,100";
	if( ! minDist ) { minDist=20; } # Arbitrary distance between reads on same "strand"
	if( ! size ) { size=2; } # Default size value - less than 2 is somewhat unreadable
	if( ! capSize ) { capSize=2; } # Size for start and end flags
	done=1;
}{
	while ( $3 !~ /chr[0-9XY]*$/ && done>0 ) { # Ignore non-canonical chr
		done=getline;
	}
	if ( $3 != chr ) { # New chromosome - reset strand count
		chr=$3;
		lastEnd=0;
		strand=1;
	}
	methCalls=substr($16, 6) # Remove leading characters in meth calls - which should always be in column 16
	calls=patsplit(methCalls, meth, /[zZ]/, len); # Split methylation calls into inter-CpG distances (len) and CpG meth calls (meth)
	if ( calls > 0 && done > 0 ) {
		start=$4;
		if ( start > lastEnd + minDist ) {  # If at least minDist away from end of last reset, start at strand 1 again
			strand=1;
		} else {
			strand++;
		}
		if( strand <=1000 ) { # Only print out 1000 strands max
			print chr, start-1, start, "start", strand, ".", 0, 0, green, capSize; # print start in green

			for ( x=1 ; x<=calls ; x++ ) {
				start+=length(len[x-1]);
				if ( meth[x] == "z" ) { # unmethylated
					colour = white;
				} else if (meth[x] == "Z" ) { # methylated
					colour = black;
				} else { # Just in case
					colour = grey;
				}
				print chr, start-1, start, meth[x], strand, ".", 0, 0, colour, size;
				start++; # Increment over Cytosine
			}
			start+=length(len[x-1]); # Add final interval
			if ( strand == 1 ) {
				lastEnd=start;
			}
			print chr, start-1, start, "end", strand, ".", 0, 0, red, capSize; # print end in red
		}
	}
}
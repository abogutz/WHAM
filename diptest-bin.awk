#! /bin/awk

# Helper script to convert Bismark BAM files methlated strand line graphs
# Made by Aaron
# Last Updated 2023-05-29

BEGIN{
	OFS="\t";
	done=1;
	if( ! thresh ) { thresh=0; }
}{
	while ( $3 !~ /chr[0-9XY]*$/ && done>0 ) { # Ignore non-canonical chr
		done=getline;
	}
	methCalls=substr($16, 6) # Get rid of leading characters
	calls=split(methCalls, len, /[zZ]/, meth); # Split methylation calls into inter-CpG distances (len) and CpG meth calls (meth)
	if ( calls > thresh && done > 0 ) {
		start=$4;
		end=start;
		for ( x=0; x<= calls; x++ ) {
			end+=length(len[x]);
			end++;
		}
		methyl=0;
		unmethyl=0;
		for ( x=1; x<= calls; x++ ) {
			if(meth[x] == "z") {
				unmethyl++;
			} else if(meth[x] =="Z") {
				methyl++;
			}
		}
		methyl=methyl/(methyl+unmethyl);
				print $3, start, end, methyl;
	}
}
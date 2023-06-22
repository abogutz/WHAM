#! /bin/awk

# Helper script to convert Bismark BAM files into bigLolly files for UCSC
# Made by Aaron
# Last Updated 2023-06-21

BEGIN{
	OFS="\t";
	chr="";
	done=1;
	read1=1;
	if( ! maxGap ) { maxGap=400; } # Treat reads with larger distances between as separate
}{
	while ( $3 !~ /chr[0-9XY]*$/ && done>0 ) { # Ignore non-canonical chr
		done=getline;
	}
	if(read1 == 1) { # First read
		storeValues();
	} else { # Second Read
		if($1 != name ) {
			print prevRead;
			storeValues();
		} else {
			methCalls2=substr($16, 6);
			if($4 < r1start) {
				if($4+length(methCalls2) > r1start) {
					methCalls2=substr(methCalls2,1,r1start-$4)
					methCalls=(methCalls2)(methCalls)
				} else {
					gap=r1start-($4+length(methCalls2));
					if(gap >= maxGap ) {
						print prevRead;
						name=$1;
						r1chr=$3;
						r1start=$4;
						methCalls=substr($16, 6); # Remove leading characters in meth calls - which should always be in column 16
						r1end=r1start+length(methCalls);
						read1=0;
						prevRead=$0;
					} else {
						stuffer="";
						for(x=1; x<=gap; x++) {
							stuffer=stuffer"+";
						}
						methCalls=(methCalls2)(stuffer)(methCalls);
					}
				}
				r1start=$4
			} else if($4< r1end) {
				methCalls2=substr(methCalls2,r1end-$4)
				methCalls=(methCalls)(methCalls2)
			} else {
				gap=$4-r1end;
				if(gap >= maxGap ) {
					print $0;
				storeValues();
			read1=0;
				} else {
					stuffer="";
					for(x=1; x<=gap; x++) {
						stuffer=stuffer"+";
					}
					methCalls=(methCalls)(stuffer)(methCalls2);
				}
			}
			$3=r1chr;
			$4=r1start;
			$16="XM:Z:"methCalls;
			print $0;
			read1=1;
		}
	}
}

function storeValues () {
	name=$1;
	r1chr=$3;
	r1start=$4;
	methCalls=substr($16, 6); # Remove leading characters in meth calls - which should always be in column 16
	r1end=r1start+length(methCalls);
	read1=0;
	prevRead=$0;
}
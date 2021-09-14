#!/usr/bin/perl

#use strict;
my($start_time,$end_time);
$start_time = localtime(time);

#This script is used to filter reads which are outside of the chromosome size limit
#Reads/peaks which are outside of the chromosome size limit are problematic for UCSC Genome Browser
#It's neccessary to filter these features before data visualization
#Way to run this script, copy script to directory of read files, run the following:
#perl Remove_Outbound_Reads.pl

#Sample input file:
# chr1	3406014	3406713	MACS_peak_1	
# chr1	3466424	3468614	MACS_peak_2	
# chr1	3504190	3505963	MACS_peak_3	
# chr1	3541097	3542564	MACS_peak_4	

#Example peak outside of chromosome size limit:
#Error message trying to make bigBed file:
#End coordinate 121843922 bigger than chr11 size of 121843856 line 5102 of G91_M1_MACS2_peaks.bed

#Sample output file:
# chr1	3406014	3406713	MACS_peak_1	
# chr1	3466424	3468614	MACS_peak_2	
# chr1	3504190	3505963	MACS_peak_3	
# chr1	3541097	3542564	MACS_peak_4	

#Do a peak count before and after to confirm filtering works.
#wc -l *.bed

#populate an array of chromosome size limits
open (file2, "../genomeIndex/mm9.chrom.sizes") or die "Could not open file2\n";
my %chr_size_limit;
my ($chromosome, $chr_size);
while (<file2>){
	my $string_2 = &Trim($_);
	#Now we have two columns to read:
	my @parts_2 = split('\t', $string_2);
	($chromosome, $chr_size) = ($parts_2[0],$parts_2[1]);
	#Populate the hash:
	$chr_size_limit{$chromosome} = $chr_size;
	}#end while loop through size limit file
	close(file2);
	#Above hash prints:
	# foreach $k (sort keys %chr_size_limit) {
		# print "$k => $chr_size_limit{$k}\n";
	# }


#Input files need to follow the file name pattern:
@id = glob("*.bed");

for my $id (@id){
	print $id."\n";
	$nline = 0;
	$id =~ s/.bed//;
	#Now $id is just the file name without extension
	open IN, $id.".bed" or die "$!";
	#open OUTbed,">".$id."_Inbound_peaks.bed" or die "$!";
	open OUTbed,">../Input_files_2/".$id."_Inbound_peaks.bed" or die "$!";
	while (<IN>){
		$nline ++;
		next if $_ !~ /^chr/ ;
		my $inline = &Trim($_);
		my @parts = split('\t', $inline);#array "parts" holds the data for each column (separated by tabs)
		my $chr = $parts[0];
 		my $start = $parts[1];
		my $end = $parts[2];
		my $name = $parts[3];

		#print $strand."\n";
		#print $chr."\n";
		#print $start."\n";
		
		
		my $line = "$chr\t$start\t$end\t$name\n";
		
		$chromosome = $chr;
		
		if ($start>0 && $start<$chr_size_limit{$chromosome} && $end<$chr_size_limit{$chromosome} && $chr !~ /random/ && $chr !~ /chrM/){
			#Want inverse match: !~ /pattern/
				#print to outfile:
				print OUTbed $line;
		}#end of if statement for error checking

	}#end of while loop within file
	print "total line\t$nline\n";
	close(IN);
    close(OUTbed);
}#end of for loop going through files in dir
		
$end_time = localtime(time);
print "Start time:$start_time\n";
print "End time:$end_time\n";

# subroutine to remove whitespace from the start and end of the string
sub Trim {
	my $string = $_[0];
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}#end of the Trim subroutine

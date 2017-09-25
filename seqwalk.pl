#!/usr/bin/perl
use warnings;
use strict;

# last modified 25 09 2017
my $initial_input_seq = 'B6_test.fasta';
my $seqs_processed = 0;


index_seq($initial_input_seq,$seqs_processed);

# align();
# find_best_hit();
# filter_alignment_score();
# convert_to_bam();

sub index_seq{
    my ($file,$seqs_processed) = @_;
    open (IN,$initial_input_seq) or die "Failed to open FastA input: $!\n";
    my $header;
    my $seq; 
    while (<IN>){
	chomp; 
	if ($_ =~ /^(>.*)/){
	    last if ($seqs_processed > 0);
	    $header = $1;
	    warn "Setting header to $1\n";
	    $seqs_processed++;
	}
	else{
	    $seq .= $_;
	}
    }
    my $outfile = "one_seq_only.fa";
    open (INDEX,">",$outfile) or die "Failed to write out FastA sequence: $!\n";
    print INDEX "$header\n$seq\n";
    close INDEX or die "Failed to close INDEX filehandle: $!\n";

    # created Last index
    warn "Now creating Last index for this sequence\n";

    system ("lastdb -v mydb $outfile");
    warn "done\n\n";
    
}


sub align {
# By default Last will report many alignments. Last-split looks for a unique best alignment for each part of each query.
# It will therefore only return the best hit for each alignment, but different parts may be aligned to different places.
# Don't use find_best_hit() if last-split is used.

	my @files = <*.fa>;

	foreach my $file (@files) {
		$file =~ s/\.fa//;
		#$file =~ s/read_files\///;
		warn "Now aligning $file \n";
		system ("lastal -P32 -r5 -q20 -a4 -b3 ../mydb $file.fa | last-split -m1e-6 > $file.maf");
	}
}





sub find_best_hit {
# Last reports many alignments as a default and this subroutine finds the one 
# with the highest alignment score. It makes no sense to use this when last-split has run.
# Also, the maf(-like) format that last-split produces is diffferent.

	my @files = <*.maf>;


	foreach my $file(@files) {



		open (my $in, $file) or die "Cannot read file: $!";
		$file =~ s/alignments\///;
		warn "$file\n";


		my $outfile = "best_hit/best_hit_$file";

		open (my $out,'>',$outfile) or die "Cannot open out file $outfile: $!";

		warn "$outfile\n";

		my $current_score;
		my $highest_score;
		my $current_ref_line;
		my $current_query_line;
		my $highest_ref_line;
		my $highest_query_line;
		my $current_alignment;
		my $highest_alignment;



		while (<$in>) {

			chomp;
			if (/^a/) {
				$current_alignment = $_;

				unless ($highest_score) {
					$highest_score = 1;
					}
				my ($score)  = (split /\s+/)[1];
				$score =~/(\d+)/;
				$current_score = $1;
				# warn "current score is $current_score\n";

				$current_ref_line = <$in>;
				$current_query_line = <$in>;
				# warn "$current_ref_line\n";
				# warn "$current_query_line\n";

				if ($current_score > $highest_score) {
					$highest_score = $current_score;

					$highest_alignment = $current_alignment;
					warn "highest score is $highest_score\n";
					$highest_ref_line = $current_ref_line;
					$highest_query_line = $current_query_line;
				}	
			}
		}
		#warn "The highest score in $file is $highest_score\n
		#$highest_alignment\n
		#$highest_ref_line\n
		#$highest_query_line\n";

		print $out "$highest_alignment\n$highest_ref_line$highest_query_line";
		close $out or die $!;
	}
}



sub filter_alignment_score {
# Although last-split and best_hit will find a unique alignment for each (part of the) query, sequences which are
# not actually present in the reference will be put in "the best wrong place". They usually have low alignments
# scores, so this subroutine can filter out alignments below a specified threshold.

	my @files = <*.maf>;

	my $threshold = 1300; # alignment score below which alignments will be excluded
	my $retained = 0;	  # number of alignments written to the filtered file
	my $kicked = 0;       # number of alignments not passing the score threshold


	
	foreach my $file(@files) {

		unless ($file =~ /_scorefiltered.maf/) {
			warn "\nFiltering $file\n";
			open (my $in, $file) or die "Cannot read file: $!";
			$file =~ s/\.maf//;
			my $outfile = $file."_scorefiltered.maf";
			open (my $out,'>',$outfile) or die "Cannot open out file $outfile: $!";
		
	
			while (<$in>) {
	
				chomp;
				if (/^#/) {
					# These lines contain general alignment info
					print $out "$_\n";
				}
	
				if (/^a\sscore/) {
					# This is the score line for each alignment
					my $score_line = $_;
					my ($score) = (split /\s+/)[1];
					$score =~/(\d+)/;
					$score = $1;
					#warn "$score_line\n";
					#warn "The score for this alignment is $score\n";
		
					if ($score > $threshold) {
			
						++ $retained;
						my $alignment = <$in>;
						chomp;
			
						my $line = <$in>;
						chomp;
						$alignment .= $line;
			
						until ($line =~ /^\n+/) {
							$line = <$in>;
							chomp;
							$alignment .= $line;
						}
						print $out "$alignment\n";
					}
				
				else {++$kicked};
				}
	
			}
			my $number_of_alignments = $retained + $kicked;
			warn "With a Last alignment score threshold of $threshold\n$retained of $number_of_alignments alignments were printed to filtered file\n$kicked of $number_of_alignments alignements were excluded\n\n\n";
	
		close $in or die $!;
		close $out or die $!;
		}
	}
}





sub convert_to_bam {

	my @files = <*_scorefiltered.maf>;
	
	foreach my $file (@files) {
		$file =~ s/\.maf//;
		warn "Now converting $file.maf to bam\n";
		system ("maf-convert sam $file.maf | samtools view -bS -t ../SAM_header.txt - > $file.bam");
		# maf-convert does not produce a SAM header that can be used to create a bam file.
		# The reference name and length can be specified in a tab delim file and accessed
		# with the -t option.
	
	
	}

}

#!/usr/bin/perl
use warnings;
use strict;

my $stepsize = 30;
my $readlength = 100;
my $number = 20;

open (IN,$ARGV[0]) or die "Please give me a .fa file!\n";

my $id = <IN>;
chomp $id;
my $seq;

while (<IN>){
    chomp;
    $seq .= $_;
}
# warn "$id\n$seq\n";

my $count = 0;
my $offset = 0;
while ($count < $number){
    ++$count;
    my $tiled_seq = substr($seq,$offset,$readlength);
    if ($count%2 ==0){
	# reverse comeplementing the sequence
	$tiled_seq = rc($tiled_seq);
    }
    print ">seq$count\n$tiled_seq\n";
    $offset += $stepsize;
}

sub rc{ # reverse complements a sequence
    my $seq = shift;
    $seq = reverse $seq;
    $seq =~ tr/CAGT/GTCA/;
    return $seq;
}

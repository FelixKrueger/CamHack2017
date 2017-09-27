#!/usr/bin/perl
use warnings;
use strict;

my $stepsize = 5000;
my $readlength = 10000;
my $number = 5;

open (IN,$ARGV[0]) or die "Please give me a .fa file!\n";

my $id = <IN>;
chomp $id;
my $seq;

while (<IN>){
    chomp;
    $seq .= $_;
}
warn "$id\n$seq\n";

my $count = 0;
my $offset = 0;
while ($count < $number){
    my $tiled_seq = substr($seq,$offset,$readlength);
    print ">seq$count\n$tiled_seq\n";
    $offset += $stepsize;
    ++$count;
}

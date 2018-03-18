#!/usr/bin/perl -w
use strict;

my($file)=@ARGV;
open(IN,"$file");

$/=">";
while(<IN>){
    chomp;   
    next if ($_ eq "");
    my @single = split("\n",$_,2);
    if ($single[0] =~ /(.*)L2/){
        print  ">".$single[0]."\n".$single[1];
    }
}

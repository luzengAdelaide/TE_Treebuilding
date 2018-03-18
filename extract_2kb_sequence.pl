#! /usr/bin/perl -w
use strict; 
my $usage="$0 <fasta> > <output>
";

my $pad = 0;  

my ($file) = @ARGV;
open(IN, $file) or die "can't open file $file:$!\n";
my $ch = '';
my $seq;
my $seq_name;
my %gs;
my $bin =0; 
my %hash;
my $head;
my $seq2;
my $head2;
my ($chr,$start,$end,$sym,$gid);
my $out = "2kb.".$file;
open(OUTA,">$out");

$/=">";
while (<IN>){
    chomp;
    next if ($_ eq "");
    $_=~s/>//;
    next if ($_ eq "");
    my @tmp = split("\n",$_,2);
    $head = $tmp[0];
    $head =~s/\>//;
    $seq2 = $tmp[1];
    $seq2 =~s/\n//g;
    print OUTA "\>$_", if length($seq2) > 2000;
}

close IN;

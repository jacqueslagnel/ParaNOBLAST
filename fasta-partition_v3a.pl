#!/usr/bin/perl -w
use strict;
use Symbol;

#modif 27/09/2014 
#jacques 21/07/2014


my $line="";
my %hit_list=();
my $total=0;
my $pos = 0;
my $toto="";
my $seq="";
my $x=0;
my @sumseqsfiles=();
my @sortedkeys=();
my %seqsl=();
my %seqs=();


#######################################################
my $num_args = $#ARGV + 1;
if ($num_args < 2) {
	print "Optimized fasta file partitioning in a given # of chunks\nmodif:07/08/2014\n";
	print "Usage: $0 <nb chunks> <fasta file>\n";
	print "output: Q_000.fas to Q_<nb chunks -1>.fas files\n";
	exit 0;
}

my $nbchunks=$ARGV[0]; # nb of chunks
my $f1f=$ARGV[1]; #in fasta
#######################################################

############## read fasta file ########################
open (SEQ,"<$f1f") or die "no file\n";
while($line=<SEQ>){
	if($line =~/^>/){
		$toto=$line;
		chomp($toto);
		$toto =~ s/^>//g;
		$seqsl{$toto}=0;
		$pos = tell();
		$seq="";
		while($line=<SEQ>){
			if ($line =~/^>/){
				seek SEQ, $pos, 0;
				last;
			}
			$pos = tell();
			chomp($line);
			$seq .=$line;
		}
		$seqs{$toto}=$seq; # hash with seq
			$seqsl{$toto}=length($seq); # hash with seq length (could be combine...)
	}
}
close (SEQ);

for ($x=0;$x<$nbchunks;$x++){
	$sumseqsfiles[$x] =0;
}


`rm -f Q_[0-9][0-9][0-9].fas`;
#---- sort hash by seq len desc ---------------------------
@sortedkeys=sort {$seqsl{$b} <=> $seqsl{$a} } keys( %seqsl);
undef %seqsl;
################## fill the files seq by seq ##############
while(@sortedkeys){
	for ($x=0;$x<$nbchunks;$x++){
		if(@sortedkeys){
			my $str=sprintf("%03d",$x);
			open(OUT, '>>', "Q_${str}.fas") or die "no out\n";
			my $kk=shift(@sortedkeys);
			print OUT ">$kk\n";
			print OUT $seqs{$kk},"\n";
			close(OUT);
			delete($seqs{$kk});
		}
	}
}

exit (0);


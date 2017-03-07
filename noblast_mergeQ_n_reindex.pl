#!/usr/bin/perl -w
use Text::ParseWords;
use Scalar::Util qw(looks_like_number);

#jac version 3 reorder sync with original fasta file

my $line="";
my $Q_old="";
my $Q="";
my @a=();
my $lid=0;
my $qid=0;
my %defs_nums=();
my $header="";
my $fastafile="";
my $num_args = $#ARGV + 1;

#-------------------------------- usage ----------------------------------------
if ($num_args != 2) {
   print "merges and reindexs blast outputs\nproduce by blast queries splitted (blast_query_splitted_pbs.pl)\nwith output option 19\n";
   print "Usage: $0 <full path of the original fasta file> <BLAST output folder (bls_outs/) containing the blast out parts (bls_outs/Q_*_bls.bls)>\n";
   exit 0;
}

$fastafile=$ARGV[0];
my $blsoutpath=$ARGV[1]; # full path of the blast outs

if(-d $blsoutpath){ #check if the dir exists
		  if(! -f "$blsoutpath/Q_000_bls.bls"){ #check if the file Q_000_bls.bls exists
					 die "ERROR No file:$blsoutpath/Q_000_bls.bls\n";
		  }
}else{
		  die "ERROR No folder:$blsoutpath\n";
}
#-------------------------------------------------------------------------------


#--------------------------- set hash with sed rank ID \t seq def -------------------
print "get seq rank ID for each def...\n";
#`grep '^>' $fastafile|grep -n '^>'|sed 's/:>/\t/g' >$blsoutpath/../defs_list_num.tsv`;
#open(TXT,"<$blsoutpath/../defs_list_num.tsv") or die "ERROR NO file num\n";
#while($line=<TXT>){
#   chomp($line);
#   if($line ne ""){
#      @a=split(/\t/,$line);
#      $defs_nums{$a[1]}=$a[0];
#   }
#}
#close(TXT);
#`rm -f $blsoutpath/../defs_list_num.tsv`;
$lid=0;
open(FAS,"<$fastafile") or die "ERROR NO file $fastafile\n";
while($line=<FAS>){
	if($line =~ /^>/){
		chomp($line);
		$line =~ s/^>//g;
		$line =~ s/\t.*//g;
	   $lid++;
		$defs_nums{$line}=$lid;
	}
}
close(FAS);

#-------------------------------------------------------------------------------

#----------------------------- merging Q chunks --------------------------------
print "Merging Q chunks ...\n";
#get header
`rm -f $blsoutpath/../tmp0.blsw`;
`head -n 1 $blsoutpath/Q_000_bls.bls >$blsoutpath/../tmp0.blsw`;
#cat all hists
#2 time [0-9]{1,}\t to avoid to get the ali for -m18 format
`grep -h -P '^[0-9]{1,}\t[0-9]{1,}\t' $blsoutpath/Q_*_bls.bls >>$blsoutpath/../tmp0.blsw`;
#-------------------------------------------------------------------------------

#-----------------------------pass 1 sync Q iD with def line ----------------------
print "Reindexing pass 1 ...\n";
open(BLS,"<$blsoutpath/../tmp0.blsw") or die "ERROR NO file $blsoutpath/../tmp0.blsw\n";
open(BLSO,">$blsoutpath/../tmp1.blsw") or die "ERROR NO file $blsoutpath/../tmp1.blsw\n";

#header
$header=<BLS>;
chomp($header);
open(HEAD,">$blsoutpath/../header.blsw") or die "ERROR NO file $blsoutpath/../header.blsw\n";
print HEAD "$header\n";
close(HEAD);

$lid=0;
while($line=<BLS>){
		  if($line =~ /^[0-9]{1,}\t[0-9]{1,}\t/){ #not to use anymore
					 chomp($line);
					 $lid++;
					 @a=split(/\t/,$line);
					 #@a=quotewords('\t', 0,$line);
					 if(scalar(@a)<32){
						die "ERROR wrong nb of fields for:$line\n";
					 }
					$a[0]=$lid;
					my $tt=$a[4];
					$tt =~ s/\"//g;
					 if(exists($defs_nums{$tt})){
						$a[1]=$defs_nums{$tt};
					}else{
						die "no Q=[$a[4]] [$tt]\n";
					}
					 $line =join("\t", @a);
					 print BLSO "$line\n";
		  }
}
close(BLS);
close(BLSO);
`rm -f $blsoutpath/../tmp0.blsw`;
#-------------------------------------------------------------------------------

#--------------------------- pass 2 sort QID,HIT,HSP ---------------------------------
print "Reindexing pass 2 ...\n";
`sort -t\$'\t' -k2,2n -k3,3n -k4,4n $blsoutpath/../tmp1.blsw >$blsoutpath/../tmp2.blsw`;
`rm -f $blsoutpath/../tmp1.blsw`;

#--------------------------- pass 3 index line ----------------------------------------
print "Reindexing pass 3 ...\n";
open(BLS,"<$blsoutpath/../tmp2.blsw") or die "ERROR NO file $blsoutpath/../tmp2.blsw\n";
open(BLSO,">$blsoutpath/../tmp3.blsw") or die "ERROR NO file $blsoutpath/../tmp3.blsw\n";

$lid=0;
while($line=<BLS>){
        if($line =~ /^[0-9]{1,}\t[0-9]{1,}\t/){ #not to use anymore
                chomp($line);
                $lid++;
                @a=split(/\t/,$line);
                if(scalar(@a)<32){
                  die "ERROR wrong nb of fields for:$line\n";
                }
               $a[0]=$lid;
                $line =join("\t", @a);
                print BLSO "$line\n";
        }
}
close(BLS);
close(BLSO);
`rm -f $blsoutpath/../tmp2.blsw`;
#-------------------------------------------------------------------------------

#----------------------------- pass 4 sort line --------------------------------
print "Reindexing pass 4 ...\n";
`sort -t\$'\t' -k1n $blsoutpath/../tmp3.blsw >$blsoutpath/../tmp4.blsw`;
`rm -f $blsoutpath/../tmp3.blsw`;
#-------------------------------------------------------------------------------

#------------------------------ add header --------------------------------------
`cat $blsoutpath/../header.blsw $blsoutpath/../tmp4.blsw >$blsoutpath/../bls_merged.bls`;
`rm -f $blsoutpath/../*.blsw`;
#-------------------------------------------------------------------------------

exit(0);


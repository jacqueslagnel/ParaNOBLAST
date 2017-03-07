#!/usr/bin/perl -w
use strict;
use Time::HiRes qw ( time alarm sleep );
use Symbol;
use IO::Handle;
use File::Basename;
$| = 1;
#############################################################
#jacques v2a
#
#added merge&indexing, bls2mysql and bls2xml in hold mode job
#12/02/2016 added extra subject def field for Blast2GO with uniprot db 
#############################################################
my $num_args = $#ARGV + 1;
if ($num_args != 8) {
	print "Run BLAST in parallel using query file splitted (version 2).\nThree BLAST outputs generated:\n1) Extended tabular (tsv) format parser free for direct blast post-analysis using for example awk (blastall -m 19 option: See NoBLast): /mypath/myblastoutput/bls_myblastoutput.bls\n2) XML format (without alignment): /mypath/myblastoutput/bls_myblastoutput.xml\n3) MySQL database (see Noblast): imbg_bls_myblastoutput (user:imbg and password:(ask me if you don't have it. Jacques)\n";
	print "\nUsage: $0 <nb jobs (queries chunks)> <fasta file full path> <out path full path> <BLAST db name full path> <blast prog [blastn|blastp|blastx]> <E-value> <nb hits to keep> <filter for low complexity: [T|F]>\n\n";
	exit 0;
}

my $jobs       = $ARGV[0];
my $file_fasta = $ARGV[1];
my $out_path   = $ARGV[2];
my $db_name    = $ARGV[3];
my $blast_prog = $ARGV[4];
my $evalue     = $ARGV[5];
my $maxhits    = $ARGV[6];
my $foutput    = $ARGV[7];
my $queue      = "batch";
#my $queue      = "fast";
my $walltime   = "60:00:00";


my $jobidlist="";
my $jobid="";
my $filename="";
my $directories="";
my $suffix="";
my $outf="";
my $seed=0;
my $tmpout="";
my $mypath="";
my $nbchunk=$jobs;
my $bls_logs="";
my $bls_jobs="";
my $bls_outs="";
my $line="";
my $seqbychunk=0;
my $counter=0;
my $x=0;
my $y=0;
my $chunkname="";
my $nbseq=0;
my $chunkpath="";
my $job="";
my @listchunk=();
my $datestr='';
my $ext='bls';
my $baseoutput="";
my $blsdb=$db_name;
my $chkdb="";

#######################################################
if ($foutput ne "F" && $foutput ne "T"){
	die "ERRO Wrong filter for low complexity option must be: F or T\n";
}

$chkdb=`ls ${db_name}.[pn]in`;
chomp($chkdb);
if( $chkdb !~ /\.[np]in/){
	die "ERROR BLAST db not found or not formated: ${db_name}\n";
}
$blsdb=~ s/.*\///g;

chomp($out_path);
$out_path=~ s/[\/]{1,}$//;
$baseoutput=$out_path;
chomp($baseoutput);

$baseoutput=~ s/.*\///g;
$baseoutput="bls_$baseoutput";
$baseoutput=~ s/\./_/g;

$chunkpath="$out_path/qchunks";
$bls_logs="$out_path/bls_logs";
$bls_jobs="$out_path/bls_jobs";
$bls_outs="$out_path/bls_outs";

if( $out_path !~ /(^~\/)|(^\/)[A-Za-z0-9._-]/){
   die "ERROR wrong path format ${out_path}\n";
}
if (length($out_path)<2){
   die "ERROR wrong path format ${out_path}\n";
}

`mkdir -p ${out_path}`;

if(! -d ${out_path}){
	die "ERROR failed to create ${out_path} folder\n";
}
if(! -W ${out_path}){
   die "ERROR the folder ${out_path} is no writable\n";
}

if(! -s $file_fasta){
	die "ERROR  fasta file not found or empty ${file_fasta}\n";
}


`rm -f ${chunkpath}/*`;
`mkdir -p ${chunkpath}`;
`rm -f ${bls_logs}/*`;
`mkdir -p ${bls_logs}`;
`rm -f ${bls_jobs}/*`;
`mkdir -p ${bls_jobs}`;
`rm -f ${bls_outs}/*`;
`mkdir -p ${bls_outs}`;


if ($queue eq "fast"){
   $walltime="03:00:00";
}


print "Running with:\nnb jobs:\t$jobs\nFasta file:\t$file_fasta\noutput path:\t$out_path\nBLAST db:\t$db_name\nBLAST program:\t$blast_prog\nE-value:\t$evalue\nnb hits 2keep:\t$maxhits\nFilter for low complexity:\t$foutput\n";
print "\nRunning on queue: \"$queue\" with max time: $walltime\n";
print "\nOutput files:\n${baseoutput}.bls\n${baseoutput}.xml\nMySQL db: imbg_${baseoutput}\n\n";


$nbseq=`grep -c '^>' $file_fasta`;
chomp($nbseq);

if($nbseq <1){
	die "ERROR no sequence found in $file_fasta\n";
}
print "Nb seq. in \"$file_fasta\": $nbseq\n";


$seqbychunk=int($nbseq / $nbchunk);
print "Nb sequence by chunk: $seqbychunk\n";

if ($seqbychunk >2){
	`rm -f Q_[0-9][0-9][0-9].fas`;
	#Optimized fasta file partitionning
	`/mnt/big/blastdb/scripts/fasta-partition_v3a.pl $jobs $file_fasta`;
	`mv Q_[0-9][0-9][0-9].fas ${chunkpath}`;
	$x=1;
	$y=0;

}else{
	print "No need to split Q fasta file\n";
	`cp $file_fasta ${chunkpath}/Q_full.fas`;
	$nbchunk=1;
}

@listchunk=`ls ${chunkpath}/Q_*.fas`;
#print "list=@listchunk\n";
$datestr='echo \`date\`';
$ext='bls';

#------------------- submit all jobs ---------------------------
$jobidlist="-W depend=afterok";
foreach my $inf (@listchunk){
	chomp($inf);
	($filename, $directories, $suffix) = fileparse($inf, qr/\.[^.]*/);
	$outf="${bls_outs}/${filename}";
	print "${filename}\n";
	$seed=int(rand(10000000))+10000000;
	$tmpout="/tmp/$filename\_$seed.bls.tmp";
$job="#!/bin/bash
#PBS -l walltime=$walltime
#PBS -q $queue
#PBS -d $out_path
#PBS -N bls_chunk_${filename}
#PBS -o ${bls_logs}/bls_${filename}.pbs.out
#PBS -j oe
#PBS -m n

${datestr}
/usr/local/bin/blastall -b ${maxhits} -v ${maxhits} -i $inf -d ${db_name} -p ${blast_prog} -F ${foutput} -e ${evalue} -m 19 -o $tmpout
${datestr}

mv -f $tmpout ${outf}_bls\.$ext
mv -f $tmpout\.par ${outf}_bls\.$ext\.par
exit 0
";

`echo "${job}" >${bls_jobs}/qsub_${filename}.sh`;
`chmod +x ${bls_jobs}/qsub_${filename}.sh`;
$jobid=`qsub ${bls_jobs}/qsub_${filename}.sh`;
chomp($jobid);
$jobid=~ s/\..*//g;
$jobidlist="$jobidlist:$jobid"
}

print "list of jobsids=$jobidlist\n";

#----------------- merge and reindex -------------------------
$job="#!/bin/bash
#PBS -l walltime=$walltime
#PBS -q $queue
#PBD -d $out_path
#PBS -N blsmerge
#PBS -o ${bls_logs}/merged_mysql_xml.out.log
#PBS -j oe
#PBS -m n

nbchunks=$jobs
echo \\\"Checking if we get all the $jobs bls chunks\\\"
for i in {1..20}
do
 nbbls=\\\`ls $bls_outs/Q_*_bls.bls.par|wc -l\\\`
 if [ \\\"\\\$nbbls\\\" -eq \\\"\\\$nbchunks\\\" ]
 then
  echo \\\"OK, got all bls chunks\\\"
  sleep 5
  break
 fi
 echo -en \\\"\\\$i,\\\"
 sleep 5
done
echo

${datestr}
echo \\\"----------- create file:${baseoutput}.bls --------------\\\"
/mnt/big/blastdb/scripts/noblast_mergeQ_n_reindex.pl ${file_fasta} ${bls_outs}
mv ${bls_outs}/../bls_merged.bls ${bls_outs}/../${baseoutput}.bls

echo \\\"----------- bls to MySQL:db =imbg_${baseoutput} --------\\\"
/mnt/big/blastdb/scripts/blastm19_to_MySQL.pl ${bls_outs}/../${baseoutput}.bls ${baseoutput}

echo \\\"----------- generates xml file:${baseoutput}.xml -------\\\"
/mnt/big/blastdb/scripts/noblast2xml_v08a.pl ${baseoutput} ${evalue} ${maxhits} ${blast_prog} ${blsdb} ${bls_outs}/../${baseoutput}.xml
#add extra field for S def when is swissprot for Blast2GO 
sed -i '/<Hit_id>/ s/<\\\/Hit_id>/| |<\\\/Hit_id>/g' ${bls_outs}/../${baseoutput}.xml
${datestr}

exit 0
";

`echo "${job}" >${bls_jobs}/qsub_noblast_mergeQ_n_reindex.pbs`;
`chmod +x ${bls_jobs}/qsub_noblast_mergeQ_n_reindex.pbs`;
$jobid=`qsub ${jobidlist} ${bls_jobs}/qsub_noblast_mergeQ_n_reindex.pbs`;
print "Final jobid=$jobid\n";

exit 0;


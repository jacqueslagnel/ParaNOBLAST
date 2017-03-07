#!/usr/bin/perl -w
use strict;
use DBI;
use Time::HiRes qw ( time alarm sleep );
use Symbol;
use IO::Handle;
use Text::ParseWords;
use Scalar::Util qw(looks_like_number);

$| = 1;
############################################################
# version 7 02/10/204
# 1) to avoid timout connection => open connection after that we got all bls parts files
# 2) added remove all \\ (with \") in db def line that could generated extra "
#######################################################
my $num_args = $#ARGV + 1;
if ($num_args < 7) {
print "Merge ,reorder and set db from the running blast (run in real time)\n";
  print "Usage: $0 <nb jobs => #BLAST db chunks)> <fasta file> <in path: raw blast outputs> <out full path> <out file base (without extention)> <db name> <nb hits to keep>\n";
  exit 0;
}

my $jobs= $ARGV[0];
my $file_fasta= $ARGV[1];
my $blast_path= $ARGV[2];
my $merged_path= $ARGV[3];
my $merged_out_base= $ARGV[4];
my $db_name = $ARGV[5];
my $maxhits=$ARGV[6];

print "Running: $0 $jobs $file_fasta $blast_path $merged_path $merged_out_base $db_name $maxhits\n";
print "Output file:\n${merged_path}/${merged_out_base}.bls\n";
#######################################################

#my $blast_path='/mnt/big/data-users/jacques/erick/blast/debug/testblast/S5_1000_blasted/raw';
#my $query_n_max=`grep -c '^>' /mnt/big/data-users/jacques/erick/blast/debug/testblast/S5_1000.fas`;
#my $merged_path="/mnt/big/data-users/jacques/erick/blast/debug/testblast";
#my $db_name = 'imbg_blastx1';
############################################################
#------------ setting 4 mysql ----------------------------
my $db_user = "imbg";
my $db_pass = 'laBos.,1968';
my $db_host = 'node0';
my $bltable = 'bl';
my $blsorted = 'blsorted';
#-----------------------------------------------------------
#------------ gloval vars ----------------------------------
my $x=0;
my $y=0;
my $query_n=1;
my $query_db ="";
my $sql_db ="";
my @erow=();
my @erow1=();
my $totalhit=0;
my $sth_db;
my $query_n_max=0;

my $filecnt=0;
my $hits=0;
my %fid2files=();
my %fid2fh=();
my @handles;
my $nbfiles=0;
my $line="";
my $indexqid=0;
my $started_time=`date`;
my $ended_time="";
my $tmpcnt=0;


`mkdir -p $merged_path`;
#-------------------- create the noblast -m 19 format output file --------------------------------------
open (TAB, ">${merged_path}/${merged_out_base}.bls") || die ("could not create file ${merged_path}/${merged_out_base}.bls\n");
print TAB "\"result ID\"\t\"Query(Q) ID\"\t\"Hit rank\"\t\"Hsp rank\"\t\"Q Name\"\t\"GI number\"\t\"SID\"\t\"Subject(S) Name\"\t\"Q size\"\t\"S size\"\t\"Alignment(Ali) size\"\t\"% Ali coverage for Q\"\t\"% Ali coverage for S\"\t\"Q Ali start\"\t\"Q Ali end\"\t\"S Ali start\"\t\"S Ali end\"\t\"Score\"\t\"BIT score\"\t\"E Value\"\t\"Identity\"\t\"% Identity\"\t\"Positive\"\t\"% Positive\"\t\"Gaps Total\"\t\"% Gaps Total\"\t\"Q Gaps\"\t\"S Gaps\"\t\"Q Frame\"\t\"S Frame\"\t\"Iteration\"\t\"Comments\"\n";

#------- get total nb seq     ---------------------------------------------
$query_n_max=`grep -c '^>' $file_fasta`;
chomp ($query_n_max);
$query_n_max =$query_n_max*1;
print "Total seq: $query_n_max\n";

#------------ for time ----------------------------------------------------
chomp ($started_time);
print "started at: $started_time\n";

#--- wait 4 all files are generated ---------------------------------------
print "waiting for bls:check nb of chunks ($jobs)...\n";
$tmpcnt=`ls $blast_path/*.bls 2>/dev/null|wc -l`;chomp ($tmpcnt);
while($tmpcnt < $jobs){
	sleep(1);
	$tmpcnt=`ls $blast_path/*.bls 2>/dev/null|wc -l`;chomp ($tmpcnt);
}
print "Ok all bls files are running...\n";
print "connect now to MySQL..\n";
#---------- connect to mysql and create db -------------------------------
my $db = "dbi:mysql:dbname=test;mysql_connect_timeout=335200;host=${db_host};mysql_socket=/var/run/mysqld/mysqld.sock";
my $c_db = DBI->connect($db, $db_user, $db_pass, { RaiseError => 1, AutoCommit => 0 }) || die "Error connecting to the database: $DBI::errstr\n";

&create_blast_db();

$sql_db = "USE  $db_name;";
$sth_db = $c_db->prepare($sql_db) or die "Couldn't prepare statement: " . $c_db->errstr;
$sth_db->execute or die "Couldn't execute statement: " . $sth_db->errstr;
$sth_db->finish;
#-------------------------------------------------------------------------


#------------------ get all bls files name ---------------------------------
@handles=`find $blast_path -nowarn -type f -name "*.bls"`;
$nbfiles=scalar(@handles);
if ($nbfiles != $jobs){die "ERROR: nb of jobs ($jobs) diff. from nb of blasts ($nbfiles)...\n";}
print "nb of blast files: $nbfiles\n";
for ($filecnt=0;$filecnt<$nbfiles;$filecnt++) {
	chomp($handles[$filecnt]);
	print "$filecnt -> $handles[$filecnt]\n";
	$fid2files{$filecnt}=$handles[$filecnt];
}
@handles=();

#--- open all blast files 4 read ------------------------------------------------
for ($filecnt=0;$filecnt<$jobs;$filecnt++){
	my $fh = gensym;
	open ($fh, "<$fid2files{$filecnt}") or die "no file: $fid2files{$filecnt} to R\n";
	push(@handles,$fh);
	$fid2fh{$filecnt}=$fh;
}

#------------------------- start to read, merge and reorder -----------------------
#---------- read all Q but not the last one and assemble --------------------------
$query_n=1;
for ($query_n=1;$query_n<$query_n_max;$query_n++){
	open TMPF,">$merged_path/qtmp.tmp" or die "failed to open $merged_path/qtmp.tmp 4 W\n";
	$filecnt=0;
	$indexqid=$query_n+1;
	$line="";
	foreach my $handle (@handles) { # read for each blast file
		while(eof($handle)){ # clear EOF
			$handle->clearerr();
			sleep (1);
			print ".";
		}
		$line=<$handle>;
		
		if(($line !~ /\n$/) && ($line ne '')){ # argh.. we dont get the full line !!!
		#while(($line !~ /\"\n$/) && ($line ne '')){
			seek($handle, -length($line), 1); #place it back
			$line="";
			print "?";
			#sleep (1);
		}
		while($line !~ /^[0-9]{1,}\t$indexqid\t/){
			if($line =~ /^[0-9]{1,}\t$query_n\t/){
				#check if we get extra \" in Q_TITLE: field 4 or S_TITLE field 7
				#print "line->|$line|<-\n";
				my $str=$line;
				chomp($str);
				#my @merdeslpit=split ('\t',$str);
				my @merdeslpit= quotewords('\t', 0,$str);
				if(scalar(@merdeslpit)!=32){
					print "ERROR: didn't get the 32 blast fields #=",scalar(@merdeslpit),"\n";
					for (my $xpoub=0;$xpoub<scalar(@merdeslpit);$xpoub++){
               if(!looks_like_number($merdeslpit[$xpoub])){
                  $merdeslpit[$xpoub] ='"'.$merdeslpit[$xpoub].'"';
               }
					print "$xpoub=[$merdeslpit[$xpoub]]\n";
            }
					die "ERROR: line ori:\n$line\n";
				}
				if( $merdeslpit[9] !~/^\d+$/ ){
					die "ERROR: S_LEN NOT a number:\n$line\n";
				}
				$merdeslpit[4]=~ s/[\"\\]//g;
				$merdeslpit[6]=~ s/[\"\\]//g;
				$merdeslpit[7]=~ s/[\"\\]//g;
				$merdeslpit[4]=~ s/\s+/ /g;
				$merdeslpit[6]=~ s/\s+/ /g;
				$merdeslpit[7]=~ s/\s+/ /g;
				for (my $xpoub=0;$xpoub<scalar(@merdeslpit);$xpoub++){
					if(!looks_like_number($merdeslpit[$xpoub])){
						$merdeslpit[$xpoub] ='"'.$merdeslpit[$xpoub].'"';
					}
				}
				#$merdeslpit[4] ='"'.$merdeslpit[4].'"';
				#$merdeslpit[6] ='"'.$merdeslpit[6].'"';
				#$merdeslpit[7] ='"'.$merdeslpit[7].'"';
				$str = join("\t", @merdeslpit );
				print TMPF "$str\n"; # fill the temp file for 1 Q
			}
			while(eof($handle)){
				$handle->clearerr();
				sleep (1);
				print ".";
			}
			$line=<$handle>;
			if(($line !~ /\n$/) && ($line ne '')){ # we dont get the full line !!!
			        seek($handle, -length($line), 1); #replace it back
			        $line="";
				#print "!";
			}
		}
		seek($handle, -length($line), 1); # place the pointer of same line back onto the filehandle
	}
	close(TMPF);
	#//////////////// do the reoder with the temp file /////////////////
	&main_reorder();
	#///////////////////////////////////////////////////////////////////
	print "Q: $query_n done\n";
}

##--------- for the last Q -------------------------------------------------------------------
#---- wait and check that the blast is finish ----------------
print "waiting for the last Q:$query_n_max\n";
#--- wait for the nb jobs are finish
$tmpcnt=`ls $blast_path/*.bls.par 2>/dev/null|wc -l`;chomp ($tmpcnt);
while($tmpcnt <$jobs){
	sleep(1);
	$tmpcnt=`ls $blast_path/*.bls.par 2>/dev/null|wc -l`;chomp ($tmpcnt);
}
#--- be sure to get the last row of the par file
$tmpcnt=`grep -l -s 'RUN_TIME_FORMATED' $blast_path/*.bls.par 2>/dev/null|wc -l`;chomp ($tmpcnt);
while($tmpcnt <$jobs){
	sleep(1);
	$tmpcnt=`grep -l -s 'RUN_TIME_FORMATED' $blast_path/*.bls.par 2>/dev/null|wc -l`;chomp ($tmpcnt);
}
sleep(1);

#----------------- doing the last one -------------------------
$query_n=$query_n_max;
$line="";
print "doing the last: $query_n\n";
open TMPF,">$merged_path/qtmp.tmp" or die "failed to open $merged_path/qtmp.tmp 4 W\n";
foreach my $handle (@handles) { # for each bls file
	while($line=<$handle>){
		if($line =~ /^[0-9]{1,}\t$query_n\t/){
				#check if we get extra \" in Q_TITLE: field 4 or S_TITLE field 7
				#print "line->|$line|<-\n";
				my $str=$line;
				chomp($str);
				#my @merdeslpit=split ('\t',$str);
				my @merdeslpit= quotewords('\t', 0,$str);
				if(scalar(@merdeslpit)!=32){
					die "ERROR: didn't get the 32 blast fields:\n$line\n";
				}
				if( $merdeslpit[9] !~/^\d+$/ ){
					die "ERROR: S_LEN NOT a number:\n$line\n";
				}
				$merdeslpit[4]=~ s/[\"\\]//g;
				$merdeslpit[6]=~ s/[\"\\]//g;
				$merdeslpit[7]=~ s/[\"\\]//g;
				#$merdeslpit[4] ='"'.$merdeslpit[4].'"';
				#$merdeslpit[6] ='"'.$merdeslpit[6].'"';
				#$merdeslpit[7] ='"'.$merdeslpit[7].'"';
				$merdeslpit[4]=~ s/\s+/ /g;
            $merdeslpit[6]=~ s/\s+/ /g;
            $merdeslpit[7]=~ s/\s+/ /g;
            for (my $xpoub=0;$xpoub<scalar(@merdeslpit);$xpoub++){
               if(!looks_like_number($merdeslpit[$xpoub])){
                  $merdeslpit[$xpoub] ='"'.$merdeslpit[$xpoub].'"';
               }
            }

				$str = join("\t", @merdeslpit );
				print TMPF "$str\n"; # fill the temp file for 1 Q
		}else{
			print "ERROR: last line: $line";
		}
	}
}
close(TMPF);
#//////////////// do the reoder with the temp file /////////////////
&main_reorder();
#///////////////////////////////////////////////////////////////////
print "assemble for the last Q:$query_n done\n";

#--------- close all files ---------------------
print "closing all files...\n";
close (TAB);
for ($filecnt=0;$filecnt<$nbfiles;$filecnt++){
	close($fid2fh{$filecnt});
}
print "normal endding..\n";

#---------- ended time ----------------------------
$ended_time=`date`;
chomp($ended_time);
print "started at: $started_time, ended at: $ended_time\n";

`rm -f $merged_path/qtmp.tmp`;
`rm -f $merged_path/temp_bls_reordered`;

exit 0;
################################### end main #######################################################

################################ subs ##############################################################

#****************************** main reorder *******************************************************
sub main_reorder(){
		#------------------------------- tmpf bls in temp db table --------------
		&create_tble_bltemp();
		#print "main_reorder\n";
		#$query_db = "TRUNCATE TABLE bltemp;";
		#$query_db = "DROP TABLE IF EXISTS bltemp;";
		#$sth_db = $c_db->prepare($query_db) or die "Couldn't prepare statement: " . $c_db->errstr;
		#$sth_db->execute or die "Couldn't execute statement: " . $sth_db->errstr;
		#$sth_db->finish;
		#&create_tble_bltemp();
		$sql_db = "LOAD DATA LOCAL INFILE  '$merged_path/qtmp.tmp' into table bltemp fields terminated by '\\t' enclosed by '\"' lines terminated by '\\n';";
		$sth_db = $c_db->prepare($sql_db) or die "Couldn't prepare statement: " . $c_db->errstr;
		$sth_db->execute or die "Couldn't execute statement: " . $sth_db->errstr;
		$sth_db->finish;
		
		#------------------------------- reorder it --------------
		open (TEMP, ">$merged_path/temp_bls_reordered") || die ("could not create file $merged_path/temp_reorder.bls\n");
		&reorderit();
		close(TEMP);
		$sql_db = "LOAD DATA LOCAL INFILE  '$merged_path/temp_bls_reordered' into table bl fields terminated by '\\t' enclosed by '\"' lines terminated by '\\n';";
		$sth_db = $c_db->prepare($sql_db) or die "Couldn't prepare statement: " . $c_db->errstr;
		$sth_db->execute or die "Couldn't execute statement: " . $sth_db->errstr;
		$sth_db->finish;
}

#******************************* reodering ************************************************************************
sub reorderit(){
    my $queryID=$query_n;
	print "Reordering query=$query_n\t";
	$x=0;
	my $hit=0;
	@erow=();
	my $nbmatch=0;

	#$query_db = "SELECT count(QUERY) FROM bltemp WHERE QUERY = $queryID AND L_S >0;";
	$query_db = "SELECT count(QUERY) FROM bltemp WHERE QUERY = $queryID AND S_TITLE not like \'%NO MATCH%\';";
	$sth_db = $c_db->prepare($query_db) or die "Couldn't prepare statement: " . $c_db->errstr;
	$sth_db->execute or die "Couldn't execute statement: " . $sth_db->errstr;
	@erow = $sth_db->fetchrow_array();
	$nbmatch=$erow[0];
	$sth_db->finish;

	$query_db = "SELECT bltemp.* FROM  bltemp WHERE QUERY = $queryID ORDER BY `E_VALUE` ASC, `BIT_SCORE` DESC, `GI` DESC;";
	#$query_db = "SELECT bltemp.* FROM  bltemp WHERE QUERY = $queryID AND S_TITLE not like \'NO MATCH\' ORDER BY `E_VALUE` ASC, `BIT_SCORE` DESC, `GI` DESC;";
	$sth_db = $c_db->prepare($query_db) or die "Couldn't prepare statement: " . $c_db->errstr;
	$sth_db->execute or die "Couldn't execute statement: " . $sth_db->errstr;

	$query_db = "SELECT bltemp.* FROM bltemp WHERE (QUERY = $queryID AND HSP > 1 AND GI like ? ) ORDER BY HSP ASC";
	my $sthdb3 = $c_db->prepare($query_db) or die "Couldn't prepare statement: " . $c_db->errstr;

	#print "all selected, Total hits=",($nbnomatch + $nbmatch),"\tmatch=$nbmatch, NO match=$nbnomatch\n";
	print "# hits=$nbmatch\n";
	while(@erow = $sth_db->fetchrow_array()) {
		if(!defined($erow[9])){print "\nERROR checkerow9=@erow\n";die "erow[9] not defined, @erow\n";} #9=s size

#1 if hit exist dont print no match
		if($nbmatch>0){
			if ($erow[3]==0){ # we dont have shit
				if ($erow[9]>0){
					$totalhit++;
					$hit++;
					$erow[0]=$totalhit;
					$erow[2]=$hit;
					&print_it();
				}#else {
				#	print "hit with Hit length ==0 and with match ??\n";
				#}
			}elsif ($erow[3]==1){ # we got shit
				$totalhit++;
				$hit++;
				$erow[0]=$totalhit;
				$erow[2]=$hit;
				&print_it();
				$sthdb3->execute($erow[5]) or die "Couldn't execute statement: " . $sthdb3->errstr;
				while(@erow1 = $sthdb3->fetchrow_array()) {
					$totalhit++;
					#$hit++; #do not inc the hit when got HSPs
					$erow1[0]=$totalhit;
					$erow1[2]=$hit;
					&print_it1();
				}
			}
		}elsif($erow[7] =~ /NO MATCH/ || $erow[9] ==0){ #print 1 only no match
			$totalhit++;
			$hit++;
			$erow[0]=$totalhit;
			$erow[2]=$hit;
			&print_it();
			last;
		}else{die "\nERROR no match but str \"NO MATCH\" not found...\n";}

		if($hit >= $maxhits){
			$hit=0;
			last;
		}
}
	$sthdb3->finish;
	$sth_db->finish;
}


#***************************** print it **************************************************************************************************************
sub print_it(){
	if (!defined($erow[30])){$erow[30]=0;}
	if (!defined($erow[31])){$erow[31]="";}
	print TAB "$erow[0]\t$erow[1]\t$erow[2]\t$erow[3]\t\"$erow[4]\"\t\"$erow[5]\"\t\"$erow[6]\"\t\"$erow[7]\"\t$erow[8]\t$erow[9]\t$erow[10]\t$erow[11]\t$erow[12]\t$erow[13]\t$erow[14]\t$erow[15]\t$erow[16]\t$erow[17]\t$erow[18]\t$erow[19]\t$erow[20]\t$erow[21]\t$erow[22]\t$erow[23]\t$erow[24]\t$erow[25]\t$erow[26]\t$erow[27]\t$erow[28]\t$erow[29]\t$erow[30]\t\"$erow[31]\"\n";
   print TEMP "$erow[0]\t$erow[1]\t$erow[2]\t$erow[3]\t\"$erow[4]\"\t\"$erow[5]\"\t\"$erow[6]\"\t\"$erow[7]\"\t$erow[8]\t$erow[9]\t$erow[10]\t$erow[11]\t$erow[12]\t$erow[13]\t$erow[14]\t$erow[15]\t$erow[16]\t$erow[17]\t$erow[18]\t$erow[19]\t$erow[20]\t$erow[21]\t$erow[22]\t$erow[23]\t$erow[24]\t$erow[25]\t$erow[26]\t$erow[27]\t$erow[28]\t$erow[29]\t$erow[30]\t\"$erow[31]\"\n";

}
sub print_it1(){
	if (!defined($erow1[30])){$erow1[30]=0;}
	if (!defined($erow1[31])){$erow1[31]="";}
	print TAB "$erow1[0]\t$erow1[1]\t$erow1[2]\t$erow1[3]\t\"$erow1[4]\"\t\"$erow1[5]\"\t\"$erow1[6]\"\t\"$erow1[7]\"\t$erow1[8]\t$erow1[9]\t$erow1[10]\t$erow1[11]\t$erow1[12]\t$erow1[13]\t$erow1[14]\t$erow1[15]\t$erow1[16]\t$erow1[17]\t$erow1[18]\t$erow1[19]\t$erow1[20]\t$erow1[21]\t$erow1[22]\t$erow1[23]\t$erow1[24]\t$erow1[25]\t$erow1[26]\t$erow1[27]\t$erow1[28]\t$erow1[29]\t$erow1[30]\t\"$erow1[31]\"\n";
	print TEMP "$erow1[0]\t$erow1[1]\t$erow1[2]\t$erow1[3]\t\"$erow1[4]\"\t\"$erow1[5]\"\t\"$erow1[6]\"\t\"$erow1[7]\"\t$erow1[8]\t$erow1[9]\t$erow1[10]\t$erow1[11]\t$erow1[12]\t$erow1[13]\t$erow1[14]\t$erow1[15]\t$erow1[16]\t$erow1[17]\t$erow1[18]\t$erow1[19]\t$erow1[20]\t$erow1[21]\t$erow1[22]\t$erow1[23]\t$erow1[24]\t$erow1[25]\t$erow1[26]\t$erow1[27]\t$erow1[28]\t$erow1[29]\t$erow1[30]\t\"$erow1[31]\"\n"

}

#********************************* sub to create blast MySQL database *********************************************
sub create_blast_db(){
	print "create db: $db_name\n";
	$sql_db = "DROP DATABASE IF EXISTS \`$db_name\`;";
	$sth_db = $c_db->prepare($sql_db) or die "Couldn't prepare statement: " . $c_db->errstr;
	$sth_db->execute or die "Couldn't execute statement: " . $sth_db->errstr;
	$sth_db->finish;

	$sql_db = "CREATE DATABASE  $db_name;";
	$sth_db = $c_db->prepare($sql_db) or die "Couldn't prepare statement: " . $c_db->errstr;
	$sth_db->execute or die "Couldn't execute statement: " . $sth_db->errstr;
	$sth_db->finish;

	$sql_db = "USE  $db_name;";
	$sth_db = $c_db->prepare($sql_db) or die "Couldn't prepare statement: " . $c_db->errstr;
	$sth_db->execute or die "Couldn't execute statement: " . $sth_db->errstr;
	$sth_db->finish;

	print "CREATE TABLE bl\n";

	$sql_db = 'CREATE TABLE  bl (
			`TOTAL_HIT` int(11) default NULL,
			`QUERY` int(11) default NULL,
			`HIT` int(11) default NULL,
			`HSP` int(11) default NULL,
			`Q_TITLE` text,
			`GI` varchar(30),
			`SID` text,
			`S_TITLE` text,
			`L_Q` int(11) default NULL,
			`L_S` int(11) default NULL,
			`SIZE_ALI` int(11) default NULL,
			`POURC_Q` double default NULL,
			`POURC_S` double default NULL,
			`Q_START` int(11) default NULL,
			`Q_END` int(11) default NULL,
			`S_START` int(11) default NULL,
			`S_END` int(11) default NULL,
			`SCORE` int(11) default NULL,
			`BIT_SCORE` double default NULL,
			`E_VALUE` double default NULL,
			`IDENTITY` int(11) default NULL,
			`POUR_IDENT` double default NULL,
			`POSITIVE` int(11) default NULL,
			`POUR_POSITIVE` double default NULL,
                        `GAPS` int(11) default NULL,
                        `POUR_GAPS` double default NULL,
                        `GAPS_Q` int(11) default NULL,
                        `GAPS_S` int(11) default NULL,
                        `SENSEQ` int(11) default NULL,
                        `SENSES` int(11) default NULL,
                        `ITERATION` int(11) default NULL,
			`COMMENTS` text default NULL,
			KEY `TOTAL_HIT` (`TOTAL_HIT`),
			KEY GI (GI(30)),
			KEY `QUERY` (`QUERY`),
			KEY `HIT` (`HIT`),
			KEY `HSP` (`HSP`)
				)
				ENGINE = MYISAM DEFAULT CHARSET=latin1;';

	$sth_db = $c_db->prepare($sql_db) or die "Couldn't prepare statement: " . $c_db->errstr;
	$sth_db->execute or die "Couldn't execute statement: " . $sth_db->errstr;
	$sth_db->finish;

	$sql_db = 'CREATE TABLE `bl_filter` (
			`TOTAL_HIT` int(11) default NULL,
			`QUERY` int(11) default NULL,
			`HIT` int(11) default NULL,
			`HSP` int(11) default NULL,
			`Q_TITLE` text,
			`GI` varchar(30),
			`SID` text,
			`S_TITLE` text,
			`L_Q` int(11) default NULL,
			`L_S` int(11) default NULL,
			`SIZE_ALI` int(11) default NULL,
			`POURC_Q` double default NULL,
			`POURC_S` double default NULL,
			`Q_START` int(11) default NULL,
			`Q_END` int(11) default NULL,
			`S_START` int(11) default NULL,
			`S_END` int(11) default NULL,
			`SCORE` int(11) default NULL,
			`BIT_SCORE` double default NULL,
			`E_VALUE` double default NULL,
			`IDENTITY` int(11) default NULL,
			`POUR_IDENT` double default NULL,
			`POSITIVE` int(11) default NULL,
			`POUR_POSITIVE` double default NULL,
                        `GAPS` int(11) default NULL,
                        `POUR_GAPS` double default NULL,
                        `GAPS_Q` int(11) default NULL,
                        `GAPS_S` int(11) default NULL,
                        `SENSEQ` int(11) default NULL,
                        `SENSES` int(11) default NULL,
                        `ITERATION` int(11) default NULL,			
			`COMMENTS` text default NULL,
			KEY `TOTAL_HIT` (`TOTAL_HIT`),
			KEY GI (GI(30)),
			KEY `QUERY` (`QUERY`),
			KEY `HIT` (`HIT`),
			KEY `HSP` (`HSP`)
				) ENGINE=MyISAM DEFAULT CHARSET=latin1;';
	$sth_db = $c_db->prepare($sql_db) or die "Couldn't prepare statement: " . $c_db->errstr;
	$sth_db->execute or die "Couldn't execute statement: " . $sth_db->errstr;
	$sth_db->finish;

	$sql_db = 'CREATE TEMPORARY TABLE `bltemp` (
			`TOTAL_HIT` int(11) default NULL,
			`QUERY` int(11) default NULL,
			`HIT` int(11) default NULL,
			`HSP` int(11) default NULL,
			`Q_TITLE` text,
			`GI` varchar(30),
			`SID` text,
			`S_TITLE` text,
			`L_Q` int(11) default NULL,
			`L_S` int(11) default NULL,
			`SIZE_ALI` int(11) default NULL,
			`POURC_Q` double default NULL,
			`POURC_S` double default NULL,
			`Q_START` int(11) default NULL,
			`Q_END` int(11) default NULL,
			`S_START` int(11) default NULL,
			`S_END` int(11) default NULL,
			`SCORE` int(11) default NULL,
			`BIT_SCORE` double default NULL,
			`E_VALUE` double default NULL,
			`IDENTITY` int(11) default NULL,
			`POUR_IDENT` double default NULL,
			`POSITIVE` int(11) default NULL,
			`POUR_POSITIVE` double default NULL,
			`GAPS` int(11) default NULL,
                        `POUR_GAPS` double default NULL,
                        `GAPS_Q` int(11) default NULL,
                        `GAPS_S` int(11) default NULL,
                        `SENSEQ` int(11) default NULL,
                        `SENSES` int(11) default NULL,
                        `ITERATION` int(11) default NULL,
			`COMMENTS` text default NULL,
			KEY `TOTAL_HIT` (`TOTAL_HIT`),
			KEY GI (GI(30)),
			KEY `QUERY` (`QUERY`),
			KEY `HIT` (`HIT`),
			KEY `HSP` (`HSP`)
				) ENGINE=MyISAM DEFAULT CHARSET=latin1;';

	$sth_db = $c_db->prepare($sql_db) or die "Couldn't prepare statement: " . $c_db->errstr;
	$sth_db->execute or die "Couldn't execute statement: " . $sth_db->errstr;
	$sth_db->finish;
}

sub create_tble_bltemp() {
	$sql_db = "DROP TABLE IF EXISTS bltemp;";
	$sth_db = $c_db->prepare($sql_db) or die "Couldn't prepare statement: " . $c_db->errstr;
	$sth_db->execute or die "Couldn't execute statement: " . $sth_db->errstr;
	$sth_db->finish;

	$sql_db = 'CREATE TEMPORARY TABLE IF NOT EXISTS `bltemp` (
			`TOTAL_HIT` int(11) default NULL,
			`QUERY` int(11) default NULL,
			`HIT` int(11) default NULL,
			`HSP` int(11) default NULL,
			`Q_TITLE` text,
			`GI` varchar(30),
			`SID` text,
			`S_TITLE` text,
			`L_Q` int(11) default NULL,
			`L_S` int(11) default NULL,
			`SIZE_ALI` int(11) default NULL,
			`POURC_Q` double default NULL,
			`POURC_S` double default NULL,
			`Q_START` int(11) default NULL,
			`Q_END` int(11) default NULL,
			`S_START` int(11) default NULL,
			`S_END` int(11) default NULL,
			`SCORE` int(11) default NULL,
			`BIT_SCORE` double default NULL,
			`E_VALUE` double default NULL,
			`IDENTITY` int(11) default NULL,
			`POUR_IDENT` double default NULL,
			`POSITIVE` int(11) default NULL,
			`POUR_POSITIVE` double default NULL,
                        `GAPS` int(11) default NULL,
                        `POUR_GAPS` double default NULL,
                        `GAPS_Q` int(11) default NULL,
                        `GAPS_S` int(11) default NULL,
                        `SENSEQ` int(11) default NULL,
                        `SENSES` int(11) default NULL,
                        `ITERATION` int(11) default NULL,			
			`COMMENTS` text default NULL,
			`ALI` mediumtext default NULL,
			KEY `TOTAL_HIT` (`TOTAL_HIT`),
			KEY GI (GI(30)),
			KEY `QUERY` (`QUERY`),
			KEY `HIT` (`HIT`),
			KEY `HSP` (`HSP`)
				) ENGINE=MyISAM DEFAULT CHARSET=latin1;';

	$sth_db = $c_db->prepare($sql_db) or die "Couldn't prepare statement: " . $c_db->errstr;
	$sth_db->execute or die "Couldn't execute statement: " . $sth_db->errstr;
	$sth_db->finish;

}



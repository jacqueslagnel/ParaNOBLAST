#!/usr/bin/perl -w
use strict;
use DBI;
use Time::HiRes qw ( time alarm sleep );
use Symbol;
use IO::Handle;
$| = 1;
############################################################

#######################################################
my $num_args = $#ARGV + 1;
if ($num_args !=2) {
print "load noblast ouptuts format -m 19 in MySQL\n";
  print "Usage: $0 [blast output m19 bls] [db name]\n";
  exit 0;
}

my $blast_outs= $ARGV[0];
my $db_name = $ARGV[1];

############################################################
#------------ setting 4 mysql ----------------------------
my $db_user = "imbg";
my $db_pass = 'laBos.,1968';
my $db_host = 'node0';
my $bltable = 'bl';

if($db_name !~ /^imbg_/){
	$db_name ="imbg_${db_name}";
}

#-----------------------------------------------------------
print "bls2MySQL database name:$db_name\n";
#------------ gloval vars ----------------------------------
my $query_db ="";
my $sql_db ="";
my $sth_db;

#---------- connect to mysql and create db -------------------------------
my $db = "dbi:mysql:dbname=test;host=${db_host};mysql_socket=/var/run/mysqld/mysqld.sock";
my $c_db = DBI->connect($db, $db_user, $db_pass, { RaiseError => 1, AutoCommit => 0 }) || die "Error connecting to the database: $DBI::errstr\n";

&create_blast_db();

$sql_db = "USE  $db_name;";
$sth_db = $c_db->prepare($sql_db) or die "Couldn't prepare statement: " . $c_db->errstr;
$sth_db->execute or die "Couldn't execute statement: " . $sth_db->errstr;
$sth_db->finish;
#-------------------------------------------------------------------------
$sql_db = "LOAD DATA LOCAL INFILE  '$blast_outs' into table bl fields terminated by '\\t' enclosed by '\"' lines terminated by '\\n' IGNORE 1 LINES;";
$sth_db = $c_db->prepare($sql_db) or die "Couldn't prepare statement: " . $c_db->errstr;
$sth_db->execute or die "Couldn't execute statement: " . $sth_db->errstr;
$sth_db->finish;

exit 0;
################################### end main #######################################################

################################ subs ##############################################################

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
TOTAL_HIT int(11) default NULL,
QUERY int(11) default NULL,
HIT int(11) default NULL,
HSP int(11) default NULL,
Q_TITLE text,
GI varchar(30),
SID text,
S_TITLE text,
L_Q int(11) default NULL,
L_S int(11) default NULL,
SIZE_ALI int(11) default NULL,
POURC_Q double default NULL,
POURC_S double default NULL,
Q_START int(11) default NULL,
Q_END int(11) default NULL,
S_START int(11) default NULL,
S_END int(11) default NULL,
SCORE int(11) default NULL,
BIT_SCORE double default NULL,
E_VALUE double default NULL,
IDENTITY int(11) default NULL,
POUR_IDENT double default NULL,
POSITIVE int(11) default NULL,
POUR_POSITIVE double default NULL,
GAPS int(11) default NULL,
POUR_GAPS double default NULL,
GAPS_Q int(11) default NULL,
GAPS_S int(11) default NULL,
SENSEQ int(11) default NULL,
SENSES int(11) default NULL,
ITERATION int(11) default NULL,
COMMENTS text default NULL,
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
TOTAL_HIT int(11) default NULL,
QUERY int(11) default NULL,
HIT int(11) default NULL,
HSP int(11) default NULL,
Q_TITLE text,
GI varchar(30),
SID text,
S_TITLE text,
L_Q int(11) default NULL,
L_S int(11) default NULL,
SIZE_ALI int(11) default NULL,
POURC_Q double default NULL,
POURC_S double default NULL,
Q_START int(11) default NULL,
Q_END int(11) default NULL,
S_START int(11) default NULL,
S_END int(11) default NULL,
SCORE int(11) default NULL,
BIT_SCORE double default NULL,
E_VALUE double default NULL,
IDENTITY int(11) default NULL,
POUR_IDENT double default NULL,
POSITIVE int(11) default NULL,
POUR_POSITIVE double default NULL,
GAPS int(11) default NULL,
POUR_GAPS double default NULL,
GAPS_Q int(11) default NULL,
GAPS_S int(11) default NULL,
SENSEQ int(11) default NULL,
SENSES int(11) default NULL,
ITERATION int(11) default NULL,
COMMENTS text default NULL,
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


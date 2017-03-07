#!/bin/bash
#jacques v8a
#modif renamed the bls output file in reorder
#added blastp for nr

################################## run parallel blast on nr/nt db splitted in parts ################################
################################# noblast v 2.1  outputs -m 19 opttion #############################################
if [ $# -lt 6 ]
then
	echo "ParaNOBlast version 6 1/11/2014: Local Optimized parallel blast against ncbi nr and nt databases"
	echo "Based on NOBlast (lagnel et al. 2009)"
	echo -en "\nUSAGE:\n"
	echo "$0 <input fasta file full path> <nr|nt> <blast program: blastx|blastp|blastn> <E-value (format e.g.: 0.01, 1e-10)> <nb hit to keep> <output full path (/mypath/myblastoutput)>"

	echo -en "\nThree BLAST outputs generated:\n1) Extended tabular (tsv) format parser free (-m 19 option See NoBLast): /mypath/myblastoutput/myblastoutput.bls\n2) XML format (without alignment): /mypath/myblastoutput/myblastoutput.xml\n3) MySQL database (see Noblast): imbg_myblastoutput (user:imbg and password:(ask me if you didn't get it. Jacques))\n"

exit 0
fi

nr_path="/mnt/big/blastdb/nr/nr"
nt_path="/mnt/big/blastdb/nt/nt"

blastallprog=/usr/local/bin/blastall
# OLD NO SID=> blastallprog=/usr/local/bin/blastall_noblast
###################################################################################################################

fasta=$1
blastdb=$2
blast_prog=$3
evalue=$4
nbhits=$5
dirout=$6

dirout=${dirout%.xml}
dirout=${dirout%/}
fout=$(basename ${dirout})

hostname -f

#------------------------- check paramaters ------------------------------------
echo "----------------------------------------------------"
echo "running with:$fasta $blastdb $blast_prog $evalue $nbhits $dirout"
echo "----------------------------------------------------"
if [ `echo $fasta|grep -c '^/'` -ne 1 ] || [ `echo $dirout|grep -c '^/'` -ne 1 ]
then
	echo "ERROR: give fasta and output files with full path..."
	exit 0
fi

if [ ! -d "$dirout" ]
then
mkdir -p $dirout
echo "Output directory created: $dirout"
else
echo "ERROE: $dirout exists"
echo "Please, delete this folder and re-run the command"
exit 0
fi

if [ ! -s $fasta ]
then
echo "ERROR fasta file no found: $fasta"
 exit 0
fi

if [ ! -d "$dirout" ]
then
echo "ERROR output directory doesnt exist: $dirout"
 exit 0
fi

cnt=`echo $fasta|grep -c -v '[A-Za-z0-9\/\.\-\_]'`
if [ $cnt -gt 0 ]
then
 echo "ERROR wrong fasta file name format..."
 exit 0
fi

cnt=`echo $dirout|grep -c -v '[A-Za-z0-9\/\.\-\_]'`
if [ $cnt -gt 0 ]
then
echo "ERROR wrong output path name format ..."
  exit 0
fi


myraw="${dirout}/raw"
bls_jobs="${dirout}/bls_jobs"
bls_logs="${dirout}/bls_logs"


#----------------- blast databases ------------------------------------
if [ "$blastdb" != "nr" ] && [ "$blastdb" != "nt" ]
then
	echo "blast database must be nr|nt"
	exit 0
fi
#--------- for nr get nb seq & db len
if [ ! -s "${nr_path}/nr.dbs" ]
then
	echo "blast nr database must be formated"
	exit 0
fi
blsdb="${nr_path}/nr.*.phr"

if [ "$blastdb" == "nr" ]
then
	if [ "$blast_prog" != "blastx" ] && [ "$blast_prog" != "blastp" ]
	then
	        echo "With nr blast program must be blastx or blastp"
      		exit 0
	fi
fi

blsdb_size=`grep '^[0-9]' ${nr_path}/nr.dbs|sed '1q;d'`
blsdb_nbseq=`grep '^[0-9]' ${nr_path}/nr.dbs|sed '2q;d'`

#-------- for nt get nb seq & db len
if [ "$blastdb" == "nt" ]
then
		if [ ! -s "${nt_path}/nt.dbs" ]
		then
       			echo "blast nt database must be formated"
       			exit 
		fi
		blsdb="${nt_path}/nt.*.nhr"
		blast_prog='blastn'
		blsdb_size=`grep '^[0-9]' ${nt_path}/nt.dbs|sed '1q;d'`
		blsdb_nbseq=`grep '^[0-9]' ${nt_path}/nt.dbs|sed '2q;d'`
fi

if [ $blsdb_size -lt 1 ] || [ $blsdb_nbseq -lt 1 ]
then
	echo "ERROR: Wrong db size |# seq"
	exit 0
fi

myjobs=`ls ${blsdb}|wc -l`
##### for reorder ----------------------------------
#myjobs=71

mydb=${fout##*/}
mydb=`echo $mydb|sed 's/\.//g'`
mydb=`echo $mydb|sed 's/ /_/g'`
mydb=`echo $mydb|sed 's/-/_/g'`
mydb="imbg_${mydb}"
echo "--------- paths results -------------------"
echo "# db chunks: $myjobs"
echo "dir out: $dirout"
echo "raw out blast: $myraw"
echo "bls_jobs: $bls_jobs"
echo "bls_logs: $bls_logs"
echo "--------- outputs -------------------------"
echo "1) BLAST output merged file -m19 fromated: ${dirout}/${fout}.bls"
echo "2) BLAST output merged xml: ${dirout}/${fout}.xml"
echo "3) BLAST MySQL db: $mydb"
echo "-------------------------------------------"

#----------------  create paths ---------------------------------------------
rm -f ${dirout}/${fout}.xml  ${dirout}/${fout}.bls
rm -f ${myraw}/bls_*
mkdir -p $myraw
outcmd=$?
if [ $outcmd -ne 0 ]
then
        echo "ERROR $outcmd: mkdir failed for: $myraw"
	exit 1
fi

rm -f ${bls_jobs}/bls_*
mkdir -p $bls_jobs
outcmd=$?
if [ $outcmd -ne 0 ]
then
        echo "ERROR $outcmd: mkdir failed for: $bls_jobs"
	exit 1
fi

rm -f ${bls_logs}/bls_*
mkdir -p $bls_logs
outcmd=$?
if [ $outcmd -ne 0 ]
then
        echo "ERROR $outcmd: mkdir failed for: $bls_logs"
	exit 1
fi


#---------------- loop over blast db parts ---------------------------------
echo "----------------- jobs -----------------------------"
for i in `ls ${blsdb}`
do

if [ "$blastdb" == "nr" ]
then
blastdbf=${i%.phr}
else
blastdbf=${i%.nhr}
fi

brut=${blastdbf##*/}
brut=`echo $brut|sed 's/\./_/g'`

mydate='echo `date`'
echo "#!/bin/bash
#PBS -l walltime=336:00:00
#PBS -N bls_${brut}
#PBS -o $bls_logs/bls_${brut}.pbs.out
#PBS -j oe
#PBS -m n

${mydate}
echo \"${blastallprog} -b ${nbhits} -v ${nbhits} -x ${blsdb_nbseq} -z ${blsdb_size} -i ${fasta} -d ${blastdbf} -p ${blast_prog} -e ${evalue} -m 19 -o $myraw/bls_${brut}.bls\"
${blastallprog} -b ${nbhits} -v ${nbhits} -x ${blsdb_nbseq} -z ${blsdb_size} -i ${fasta} -d ${blastdbf} -p ${blast_prog} -e ${evalue} -m 19 -o $myraw/bls_${brut}.bls
${mydate}

exit 0
" >${bls_jobs}/bls_job_${brut}.sh

done
#----------------------- end loop ------------------------------------------

echo "sending to the queue...."

#----------------------- in realtime run the reodering ---------------------

echo "#!/bin/bash
#PBS -l walltime=500:00:00
# PBS -q bigmem
#PBS -N bls_Reorder
#PBS -o $bls_logs/bls_reorder.pbs.out
#PBS -j oe
#PBS -m n
# PBS -l nodes=node1:ppn=1
#PBS -l nodes=1:ppn=1

echo \"----------------------- in realtime run the reodering ---------------------\"
${mydate}
echo \"/mnt/big/blastdb/scripts/blast_reorder_v8a.pl ${myjobs} ${fasta} ${myraw} ${dirout} ${fout} ${mydb} ${nbhits}\"
/mnt/big/blastdb/scripts/blast_reorder_v8a.pl ${myjobs} ${fasta} ${myraw} ${dirout} ${fout} ${mydb} ${nbhits}
echo \"blast and reordering ended .....at:\"
${mydate}

echo \"----------------- when finished start xml file generation from blast db -----------------------------------\"
echo \"starting noblast to 1 file xml\"
echo \"/mnt/big/blastdb/scripts/noblast2xml_v08a.pl ${mydb} ${evalue} ${nbhits} $blast_prog $blastdb ${dirout}/${fout}.xml\"
/mnt/big/blastdb/scripts/noblast2xml_v08a.pl ${mydb} ${evalue} ${nbhits} $blast_prog $blastdb ${dirout}/${fout}.xml

echo \"------------ ended ---------------------------------\"
${mydate}

exit 0
" > ${bls_jobs}/bls_Reorder_xml_job.sh

#----------------------- submit blast jobs to the queue --------------------
chmod +x ${bls_jobs}/*.sh
echo ""
#echo "----- First the reordering in background ---------"
#echo "qsub ${bls_jobs}/bls_Reorder_xml_job.sh"
qsub ${bls_jobs}/bls_Reorder_xml_job.sh >/dev/null
#echo ""
#echo "------ and after all blast jobs ------------------"
for i in `ls ${bls_jobs}/bls_job_*.sh`
do
	#echo "qsub job: $i"
	qsub ${i} >/dev/null
done
echo "OK: sent ${myjobs}+1 jobs to the batch queue...."

exit 0


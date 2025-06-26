##$1=species GS1_002
##$2=assemblies, e.g. /data9/home/zlgu/project/NCII/assembly/Nip.fa(IRGSP-1.0_genome.fasta)
##$3=prefix of RNA-seq Data, e.g. /data9/home/zlgu/project/NCII/assembly/annotation/rawdata/Nip
##运行前注意核对基因组文件和RNA-seq文件的命名规则
#################################################
##genome=*.fa (assembly/)
##RNA-seq=*FRINGE_clean_R1.fastq,*FRINGE_clean_R2.fastq (annotation/rawdata)
##*LEAF_clean_R1.fastq,*LEAF_clean_R2.fastq
##*ROOT_clean_R1.fastq,*ROOT_clean_R2.fastq
##*SEEDLING_clean_R1.fastq,*SEEDLING_clean_R2.fastq
#################################################
##environment: genometool
#################################################
##command line:sh anno_step1.sh GS1_002 /data9/home/zlgu/project/NCII/assembly/GS1-002_p.filter.fa /data9/home/zlgu/project/NCII/assembly/annotation/rawdata/final_trim/1GS-002
#################################################

######RepeatMasker######
mkdir /data9/home/zlgu/project/NCII/assembly/annotation/$1
proc="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$1".processing.log"
echo "Run RepeatMasker (hard)..." >$proc
/data6/tool/RepeatMasker_4.0.6/RepeatMasker -engine rmblast -pa 20 -nolow -species rice -dir /data9/home/zlgu/project/NCII/assembly/annotation/$1 $2
temp1=$2
temp2=${temp1%.fa}
temp3=${temp2##*/}
temp4="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$temp3".fa.masked"
temp5="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$temp3".masked.fa"
mv $temp4 $temp5
tbl1="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$temp3".fa.tbl"
tbl2="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$temp3".masked.fa.tbl"
mv $tbl1 $tbl2
echo $temp3 >>$proc

#echo "Run RepeatMasker (soft)..." >> $proc
#/data9/home/zlgu/tool/RepeatMasker-master/RepeatMasker -e rmblast -pa 64 -nolow -species rice -xsmall -dir /data9/home/zlgu/project/NCII/assembly/annotation/$1 $2
#temp6="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$temp3".softmasked.fa"
#tbl3="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$temp3".softmasked.fa.tbl"
#mv $temp4 $temp6
#mv $tbl1 $tbl3
######Hisat2 index######
temp7=${temp5%.fa}
hisat2-build $temp5 $temp7
#######Hisat2 alignment#######
temp9=$3"FRINGE_clean_R1.fastq"
temp10=$3"FRINGE_clean_R2.fastq"
temp11="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$1"FRINGE.sam"
temp12="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$1"FRINGE.bam"
temp13="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$1"FRINGE_sort.bam"
echo "Run Hisat2 (FRINGE)..." >>$proc
hisat2 -x $temp7 -1 $temp9 -2 $temp10 --dta --rg-id fringe --rg "SM:fringe,LB:lib1,PL:ILLUMINA" -S $temp11 -p 32 && samtools view -Sb $temp11 >$temp12 && rm $temp11 && samtools sort -@ 10 $temp12 -o $temp13 && rm $temp12

temp14=$3"LEAF_clean_R1.fastq"
temp15=$3"LEAF_clean_R2.fastq"
temp16="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$1"LEAF.sam"
temp17="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$1"LEAF.bam"
temp18="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$1"LEAF_sort.bam"
echo "Run Hisat2 (LEAF)..." >>$proc
hisat2 -x $temp7 -1 $temp14 -2 $temp15 --dta --rg-id leaf --rg "SM:leaf,LB:lib1,PL:ILLUMINA" -S $temp16 -p 32 && samtools view -Sb $temp16 >$temp17 && rm $temp16 && samtools sort -@ 10 $temp17 -o $temp18 && rm $temp17

temp19=$3"ROOT_clean_R1.fastq"
temp20=$3"ROOT_clean_R2.fastq"
temp21="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$1"ROOT.sam"
temp22="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$1"ROOT.bam"
temp23="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$1"ROOT_sort.bam"
echo "Run Hisat2 (ROOT)..." >>$proc
hisat2 -x $temp7 -1 $temp19 -2 $temp20 --dta --rg-id root --rg "SM:root,LB:lib1,PL:ILLUMINA" -S $temp21 -p 32 && samtools view -Sb $temp21 >$temp22 && rm $temp21 && samtools sort -@ 10 $temp22 -o $temp23 && rm $temp22

temp24=$3"SEEDLING_clean_R1.fastq"
temp25=$3"SEEDLING_clean_R2.fastq"
temp26="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$1"SEEDLING.sam"
temp27="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$1"SEEDLING.bam"
temp28="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$1"SEEDLING_sort.bam"
echo "Run Hisat2 (SEEDLING)..." >>$proc
hisat2 -x $temp7 -1 $temp24 -2 $temp25 --dta --rg-id seedling --rg "SM:seedling,LB:lib1,PL:ILLUMINA" -S $temp26 -p 32 && samtools view -Sb $temp26 >$temp27 && rm $temp26 && samtools sort -@ 10 $temp27 -o $temp28 && rm $temp27
######samtools merge######
temp29="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$1"Merge.bam"
echo "Run merge..." >>$proc
samtools merge -@ 32 -o $temp29 $temp13 $temp18 $temp23 $temp28
######Stringtie######
temp30="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$1"Merge.gtf"
echo "Run stringtie..." >>$proc
stringtie -o $temp30 -p 32 $temp29

##merged gtf move to PASA path
temp31="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/PASA"
mkdir $temp31
temp32="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/PASA/"$1"Merge.gtf"
mv $temp30 $temp32

temp33="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/PASA/"$temp3".masked.fa"
cp $temp5 $temp33

temp34="/data9/home/zlgu/project/NCII/assembly/annotation/PASA/pasa.alignAssembly.txt"
temp35="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/PASA/pasa.alignAssembly.txt"
cp $temp34 $temp35

status="/data9/home/zlgu/project/NCII/assembly/annotation/"$1"/"$1".status.log"
echo "Step1 are Done ..." >$status

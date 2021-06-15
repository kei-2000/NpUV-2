#!/bin/bash

ref=/home/matsuo/Desktop/test/ref/GCF_000026045.1_ASM2604v1_genomic.fna
wd=/home/matsuo/Desktop/test/210614/coverage5

fastqd=/home/matsuo/Desktop/test/210609/fasta_20_link #fastqファイルの場所

F=0.3
C=5 #bwaのオプション

fastq=$fastqd/*1.fq.gz

mkdir $wd/ref
cp $ref $wd/ref/reference.fna
f=$wd/ref/reference.fna

bwa index $f

mkdir $wd/bam
mkdir $wd/vcf
mkdir $wd/fastq_QC
mkdir $wd/vcf.gz

ls $fastq | while read line

do

r1=${line%1.fq.gz}1.fq.gz
r2=${line%1.fq.gz}2.fq.gz

a=${r1#$fastqd/}
b=${r2#$fastqd/}

fastp -i $r1 \
-I $r2 \
-o $wd/fastq_QC/${a%.fq.gz}_QC.fq \
-O $wd/fastq_QC/${b%.fq.gz}_QC.fq \
-f 1 -F 1 -t 1 -T 1 -w 16 -q 30

R1=$wd/fastq_QC/${a%.fq.gz}_QC.fq
R2=$wd/fastq_QC/${b%.fq.gz}_QC.fq

BAM=$wd/bam/${a%1.fq.gz}.bam
VCF=$wd/vcf/${a%1.fq.gz}.vcf
#VCFGZ=$wd/vcf.gz/${a%1.fq.gz}.vcf

bwa mem \
-R "@RG\tID:${a%.fq.gz}\tLB:library\tSM:${a%.fq.gz}\tPL:ILLUMINA" \
-t 20 \
$f $R1 $R2 \
|elprep filter /dev/stdin $BAM \
--output-type bam --filter-unmapped-reads --filter-mapping-quality 30 --mark-duplicates --remove-duplicates --sorting-order coordinate &&\
samtools index $BAM
freebayes -b $BAM \
-f $f -F $F -C $C \
>$VCF &&\

cp $VCF $wd/vcf.gz/ & #&&\

#bgzip $VCFGZ &&\
#tabix -p vcf $VCFGZ.gz &

done

#wd=/home/kei/デスクトップ/210609/test #作業ディレクトリ指定
#vcf=/home/kei/デスクトップ/210609/test  #vcfの場所

i=0
k=1
l=1
m=0 #変数定義

#######################################vcfのマージ################################################
vmerge="bcftools merge --merge all --force-samples"

#mkdir $wd/vcf.gz
mkdir $wd/mergedvcf
mkdir $wd/vaf
#cp $vcf/*.vcf $wd/vcf.gz
ls $wd/vcf.gz/*.vcf >$wd/tmp1.txt

echo "vcfファイルの圧縮開始"

while read line;do
bgzip $line
tabix -p vcf $line.gz
echo "圧縮... $line → $line.gz"
vmerge="$vmerge $line.gz"

l="$l,$((17+$i*7)),$((20+$i*7))"
i=$(($i+1))

done <$wd/tmp1.txt

echo 完了
echo 合成vcfファイルの作成を開始

$vmerge >$wd/mergedvcf/merged.vcf

echo 完了

#######################################vcfのマージ################################################

#######################################variant allele frequencyの計算################################################

echo VAFの計算を開始

sed -e s/:/"\t"/g $wd/mergedvcf/merged.vcf >$wd/mergedvcf/sed.txt
grep -v "#" $wd/mergedvcf/sed.txt |\
cut -f $l >$wd/mergedvcf/cut.txt
cat $wd/mergedvcf/cut.txt | while read line;do

j=2

while [ $j -le $(($i*2)) ];do

a=$(echo $line |awk "la=$(($j+1))"'{print $la}')
b=$(echo $line |awk "lb=$j"'{print $lb}')

if [[ $a == *,* ]];then

list1=(${a//,/" "})
n=$(echo ${#list1[*]})

while [ $m -lt $n ];do
if [[ ${list1[$m]} == . ]];then

list4="$list4"s"0%" 
m=$(($m+1))

else

y=$((${list1[$m]}*100/$b))
list4="$list4"s""$y"%"
m=$(($m+1))

fi
done

list3="$list3"t"($list4)"
list4=""
m=0

else
if [[ $a == . ]];then

list3="$list3"t"0%" 

else

y=$(($a*100/$b))
list3="$list3"t""$y"%"

fi
fi

j=$(($j+2))

done

echo ${list3[*]} >>$wd/mergedvcf/tmp2.txt
list3=""
k=$(($k+1))

done

sed -ie s/t/"\t"/g $wd/mergedvcf/tmp2.txt
sed -ie s/%s/%,/g $wd/mergedvcf/tmp2.txt
sed -ie s/s//g $wd/mergedvcf/tmp2.txt

grep -v "#" $wd/mergedvcf/merged.vcf >$wd/mergedvcf/tmp3.txt
sed -ie s/"\t"/!/g $wd/mergedvcf/tmp3.txt
cut -d ! -f 1,2,4-6 $wd/mergedvcf/tmp3.txt >$wd/mergedvcf/left.txt
sed -ie s/!/"\t"/g $wd/mergedvcf/left.txt

grep -v "##" $wd/mergedvcf/merged.vcf >$wd/mergedvcf/tmp4.txt
sed -ie s/"\t"/!/g $wd/mergedvcf/tmp4.txt
head -n 1 $wd/mergedvcf/tmp4.txt >$wd/mergedvcf/tmp5.txt
cut -d ! -f 1,2,4-6,8,10- $wd/mergedvcf/tmp5.txt >$wd/mergedvcf/head.txt
sed -ie s/!/"\t"/g $wd/mergedvcf/head.txt

paste $wd/mergedvcf/left.txt $wd/mergedvcf/tmp2.txt >$wd/mergedvcf/tmp6.txt
cat $wd/mergedvcf/head.txt $wd/mergedvcf/tmp6.txt >$wd/vaf/vaf.tsv

echo 完了

#######################################variant allele frequencyの計算################################################
echo 一時ファイルを削除

rm $wd/tmp1.txt
rm $wd/mergedvcf/sed.txt
rm $wd/mergedvcf/cut.txt
rm $wd/mergedvcf/tmp*
rm $wd/mergedvcf/left.txt
rm $wd/mergedvcf/head.txt

echo 完了

exit 0


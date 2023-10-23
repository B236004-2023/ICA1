#!/bin/bash

#extracts columns 6 to 7 from the file 'Tco2.fqfiles' and writes the output to a file named 'srr.list'
cut -f6-7 Tco2.fqfiles >srr.list
#resd each line from 'srr.list'assigning the values in columns 6 and 7 to the variables 'old' and 'new'
#and add the fastqc command into 'Para.sh'file
while read old new
do echo "fastqc -t 4 -o ../fastqc $old $new" >Para.sh
done<srr.list

#change the format
sed -i "s/\.gz//g" srr.list
sed -i "1d"  srr.list
head srr.list
#
mkdir index
#build genome index to align the read pairs
bowtie2-build TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz index/aaa
#
cd fastq
#read lines of ssr.list and compare -1 and -2 files and output them into align.sh with sam format
while read old new
do echo "bowtie2 -x ../index/aaa -1 $old -2  $new -S ${old/_1.fq/}.sam"
done>align.sh<srr.list
#check the align.sh
head align.sh
#run loops
sh align.sh
#change the format from sam to bam
for i in *.sam; do samtools sort -o ${i/sam/}bam -O bam $i;done


#copy the bedfile and change the format from bam to bed
cp /localdisk/data/BPSM/ICA1/TriTrypDB46_TcongolenseIL3000_2019.bed ../
for i in *.bam; do bedtools bamtobed -i $i >${i/bam/}bed;done
#use for loops to reports the overlapping intervals and redirected to new files end with '.intersect.bed'
for i in *.bed
do bedtools intersect -wa -wb -a ../TriTrypDB-46_TcongolenseIL3000_2019.bed -b $i >${i/.bed/}.intersect.bed
done
#use for loops to count the intersect occurrences and output the results to '.count.txt'files
for b in *.intersect.bed
do for i in `cut -f4 $b|sort|uniq`
do echo -e "$i\t\c"
grep -w "$i" $b|wc -l
done >${b/.bed/}.counts.txt
done

#merge the table
nano average.awk
#!/usr/bin/awk -f

BEGIN {
    OFS = "\t"
    header = "Gene"
    num_files = ARGC - 1
}

{
    header = header OFS ARGV[ARGIND]
    while (getline < ARGV[ARGIND]) {
        total[$1][ARGIND] = $2
    }
    close(ARGV[ARGIND])
}

END {
    print header
    for (key in total) {
        for (i = 1; i <= num_files; i++) {
            if (!total[key][i]) {
                total[key][i] = 0
            }
        }
        printf "%s", key
        for (i = 1; i <= num_files; i++) {
            printf "%s%s", OFS, total[key][i]
        }
        print ""
    }
}

#according to the treatment induce\sample_ type\time
for C in Clone1 Clone2 WT
do for t in 0 24 48
do for T in Induced Uninduced
do echo "grep -w \"$C\" Tco2.fqfiles |grep -w \"$t\"|grep -w \"$T\"|cut -f1>$C-$t-$T.txt"
done;done;done>>all.txt
#
sh all.txt#
wc -l *nduced.txt
#
for i in *nduced.txt
do command="awk -f average.awk "
 for j in `cat $i`
 do command+=" $j.inetrsect.counts.txt "
 done ;echo "$command"
done

#loops of average
compare_1.awk
{
    sum = 0
    count = 0
    for (i = 2; i <= NF; i++) {
        sum += $i
        count++
    }
    if (count > 0) {
        average = sum / count
    } else {
        average = 0
    }
    printf "%s\t%.2f\n", $1, average
}

#calculate the average
for i in *.merge.txt;do awk -f    compare_1.awk $i >${i/.txt/}.averge.txt;done
awk -f average.awk *.merge.averge.txt

#except the first column, every two columns compare each other
compare.awk
BEGIN {FS = "\t";OFS = "\t"}

NR == 1 {
    for (i = 2; i <= NF; i++) {
        sample[i - 1] = $i
    }
    # print
    printf "Gene"
    for (i = 1; i < NF; i++) {
        for (j = i + 1; j <= NF; j++) {
            printf "\t%s vs %s", sample[i], sample[j]
        }
    }
    print ""
}

NR > 1 {
    gene = $1
    printf "%s", gene
    for (i = 2; i < NF; i++) {
        count1 = $i
        for (j = i + 1; j <= NF; j++) {
            count2 = $j
            if (count1 != 0 && count2 != 0) {
                division_result = count1 / count2
                printf "\t%.2f", division_result
            } else {
                printf "\tNA"
            }
        }
    }
    print ""
}

#compare the counts,'2'is an example of sorting
awk -f compare.awk all.average.txt|sort -nk2>group_wise.txt


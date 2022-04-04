# filter reditools output

# Written by Enrico Mossoto
# Ported by Jack Humphrey
# 2021

INPUT=$1
OUTPUT=$2
#set a minimum depth for the locus of 5
zless $INPUT | 
awk 'BEGIN{OFS=FS="\t"}{if($5>5) print $0}' |
#split multi editing events in multiple lines
awk 'BEGIN{OFS=FS="\t"} {if(length($8)>2) {split($8,a," "); for (i in a){$8=a[i];print $0}} else print $0}' |
#Format for annovar and reverse ref/alt for reverse transcripts
awk 'BEGIN{OFS=FS="\t"} {if(NR==1) true; else if (NR>1){ split($8,a,""); if ($4==0) {if ($3=="A") {$3="T"} else if ($3=="C") {$3="G"} else if ($3=="G") {$3="C"} else if ($3=="T") {$3="A"}; if (a[2]=="A") {a[2]="T"} else if (a[2]=="C") {a[2]="G"} else if (a[2]=="G") {a[2]="C"} else if (a[2]=="T") {a[2]="A"}}; print $1,$2,$2,$3,a[2],$4,$5,$6,$7,$8}}' > $OUTPUT



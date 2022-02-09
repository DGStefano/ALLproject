#!/bin/bash

#Scritp Variables :
#                   1. $1 data's path
#                   2. $2 output path
#
#                   if -v is passed script skips Sharing Index calc
#

SHARING=1

while [ -n "$1" ]; do # while loop starts

    case "$1" in

    -v) SHARING=0 ;;

    --)
        shift # The double dash makes them parameters

        break
        ;;

    *) echo "Option $1 not recognized" ;;

    esac

    shift

done

#--------------------------------------------------------------------------------------------------------------------------------
#|                                                                                                                              |
#|                                                            DEDULICAZIONE DELLE READS                                         |
#|                                                                                                                              |
#--------------------------------------------------------------------------------------------------------------------------------

echo "-------------------------------------------------------------------------------------------------------------------------------"
echo "|                                                                                                                             |"
echo "|                                                     READS DEDUPLICATION                                                     |"
echo "|                                                                                                                             |"
echo "-------------------------------------------------------------------------------------------------------------------------------"



#create output folder for dedu files
mkdir "$2/dedu_files"




#definition of paths for deduplication
tempdir="/tmp/"
wdir="$1/*.bam"
out_dedu="$2/dedu_files"



for complete_path in $wdir; do
    filename=${complete_path##*/}
    filename_noext=${filename%.*}
    
    source /home/stefano/miniconda3/etc/profile.d/conda.sh
    conda activate gatk

    gatk MarkDuplicates --TMP_DIR $tempdir --REMOVE_DUPLICATES true -I $complete_path -O $out_dedu/dedu_$filename_noext.bam -M $out_dedu/dedu_$filename_noext.txt

    source /home/stefano/miniconda3/etc/profile.d/conda.sh
    conda activate samtools
    
    samtools sort -T $tempdir/$filename_noext.sorted -o $out_dedu/sorted_dedu_$filename_noext.bam $out_dedu/dedu_$filename_noext.bam
    echo "indexing"
    samtools index $out_dedu/sorted_dedu_$filename_noext.bam

    rm "$out_dedu/dedu_$filename_noext.bam" "$out_dedu/dedu_$filename_noext.txt"
done




echo;echo;echo;echo;echo;echo

#--------------------------------------------------------------------------------------------------------------------------------
#|                                                                                                                              |
#|                                                     CHIAMATA DEI PICCHI                                                      |
#|                                                                                                                              |
#--------------------------------------------------------------------------------------------------------------------------------

#create output folder for macs2
mkdir "$2/macs2"


echo "-------------------------------------------------------------------------------------------------------------------------------"
echo "|                                                                                                                             |"
echo "|                                                        PEAKS CALLING                                                        |"
echo "|                                                                                                                             |"
echo "-------------------------------------------------------------------------------------------------------------------------------"




#definition of paths macs and macs output
tempdir="/tmp/"
wdir="$2/dedu_files/*.bam"
out_macs="$2/macs2"

for complete_path in $wdir; do
    filename=${complete_path##*/}
    filename_noext=${filename%.bam}
    
    source /home/stefano/miniconda3/etc/profile.d/conda.sh
    conda activate ifo-base
    
    #Chiamata con macs2
    macs2 callpeak -t $complete_path --tempdir $tempdir --format AUTO --nomodel --shift -100 --extsize 200 -B --SPMR --call-summit -q 0.01 -g hs --outdir $out_macs -n $filename_noext
    sort -k1,1 -k2,2n $out_macs/"$filename_noext"_treat_pileup.bdg > $out_macs/"$filename_noext"_treat_pileup_S.bdg
    
    source /home/stefano/miniconda3/etc/profile.d/conda.sh
    conda activate bedtobw
    
    #converto bed a bigwig    
    bedGraphToBigWig $out_macs/"$filename_noext"_treat_pileup_S.bdg $2/genome_info/hg19.chrom.sizes $out_macs/"$filename_noext"_treat_pileup_S.bw
    rm $out_macs/"$filename_noext"_treat_pileup_S.bdg
    rm $out_macs/"$filename_noext"_treat_pileup.bdg

    echo $filename_noext "done"
done



echo;echo;echo;echo;echo;echo
#--------------------------------------------------------------------------------------------------------------------------------
#|                                                                                                                              |
#|                                           ELIMINAZIONE DELLE BLACKLIST REGIONS                                               |
#|                                                                                                                              |
#--------------------------------------------------------------------------------------------------------------------------------

mkdir "$2/filtered_peak"


echo "-------------------------------------------------------------------------------------------------------------------------------"
echo "|                                                                                                                             |"
echo "|                                                 REMOVING BLACKLIST REGIONS                                                  |"
echo "|                                                                                                                             |"
echo "-------------------------------------------------------------------------------------------------------------------------------"


#change blacklist path for different project
blacklist="/home/stefano/Documents/IFO/folgiero/blacklistregions_ENCODE/consensusBlacklist_hg19.bed" 

narrow_peak="$2/macs2/*.narrowPeak"
out_filtered="$2/filtered_peak"

for complete_path in $narrow_peak ; do
    filename=${complete_path##*/}				
    filename_noext=${filename%.*.*}
    echo "## filtering $filename_noext ##"
    
    source /home/stefano/miniconda3/etc/profile.d/conda.sh
    conda activate ifo-base
    
    bedtools intersect -a $complete_path -b $blacklist -v > $out_filtered/filtered_$filename_noext.bed



    echo "## $filename filtered! ##"
    echo "#"
    echo "#"
    echo "#"
    echo "#"
    echo "#"
done


echo;echo;echo;echo;echo;echo

#--------------------------------------------------------------------------------------------------------------------------------
#|                                                                                                                              |
#|                                            CONTA DELLE READS PER OGNI PICCO                                                  |
#|                                                                                                                              |
#--------------------------------------------------------------------------------------------------------------------------------

mkdir "$2/multicovs"



echo "-------------------------------------------------------------------------------------------------------------------------------"
echo "|                                                                                                                             |"
echo "|                                                   READS COUNT PER PEAKS                                                     |"
echo "|                                                                                                                             |"
echo "-------------------------------------------------------------------------------------------------------------------------------"

#path for read's count
filtered_files="$2/filtered_peak/*.bed"
dedu_files="$2/dedu_files/sorted_*.bam"
out_bed="$2/multicovs"



#cicling on filtered peaks and deduplicated files
for complete_path_dedu in $dedu_files; do
    filename_compl_dedu=${complete_path_dedu##*/}	
    filename_compl_dedu2=${filename_compl_dedu##sorted_dedu_}	
    filename_noext_dedu=${filename_compl_dedu2%.*.*}

    #cicling on filtered peacks
    for complete_path_filtered in $filtered_files; do			
        filename_compl_filtered=${complete_path_filtered##*/}	
        filename_compl_filtered2=${filename_compl_filtered##filtered_sorted_dedu_}	
        filename_noext_filtered=${filename_compl_filtered2%.*}

        #matching same sample files on filtered3 and dedu3 
        if [ "$filename_noext_dedu" == "$filename_noext_filtered" ] ; then
            echo "bed su $filename_noext_dedu --- $filename_noext_filtered" ; echo
            
            
            source /home/stefano/miniconda3/etc/profile.d/conda.sh
            conda activate ifo-base
            
            #bedtools multicov between files from same sample
            bedtools multicov -bams $complete_path_dedu -bed  $complete_path_filtered > $out_bed/multicovs_filtered_sorted_dedu_"$filename_noext_filtered".fastq.bed
            echo $"done" ; echo ; echo
        fi
    done
done


echo;echo;echo;echo;echo;echo

#--------------------------------------------------------------------------------------------------------------------------------
#|                                                                                                                              |
#|                                            AGGIUNTA COLONNA CON LUNGHEZZA DEI PICCHI                                         |
#|                                                                                                                              |
#--------------------------------------------------------------------------------------------------------------------------------

mkdir "$2/multicovs_length"


echo "-------------------------------------------------------------------------------------------------------------------------------"
echo "|                                                                                                                             |"
echo "|                                                   ADDING PEAKS' LENGTH                                                      |"
echo "|                                                                                                                             |"
echo "-------------------------------------------------------------------------------------------------------------------------------"


multicovs_file="$2/multicovs/*.bed"
out_bed="$2/multicovs_length/"

script_path="/home/stefano/Documents/IFO/folgiero/script/sirianalysis/calc_peack_length.R"

for complete_path in $multicovs_file ; do
    filename_compl=${complete_path##*/}

    echo "$complete_path"
    echo "calc lenght on $filename_compl"; echo
    Rscript $script_path $complete_path $out_bed
    echo "done"
done


echo;echo;echo;echo;echo;echo

#--------------------------------------------------------------------------------------------------------------------------------
#|                                                                                                                              |
#|                                            NORMALIZZAZIONE DEI PICCHI                                                        |
#|                                                                                                                              |
#--------------------------------------------------------------------------------------------------------------------------------

mkdir "$2/nscore_peak"


echo "-------------------------------------------------------------------------------------------------------------------------------"
echo "|                                                                                                                             |"
echo "|                                                   PEAKS NORMALIZATION                                                       |"
echo "|                                                                                                                             |"
echo "-------------------------------------------------------------------------------------------------------------------------------"




# Nscore= ((peak read count / peak size)⋅106))* 103 /total mapped reads
# uso script R normalize_peak.R che normalizza e sorta

multicovs_length_file="$2/multicovs_length/*.bed"
out_bed="$2/nscore_peak/"
script_path="/home/stefano/Documents/IFO/folgiero/script/sirianalysis/normalize_peak.R"
dedu_file="$2/dedu_files/sorted_*.bam"

for complete_path_multicovs in $multicovs_length_file; do
    filename_compl_multicov=${complete_path_multicovs##*/}
    filename_compl_multicov2=${filename_compl_multicov##length_multicovs_filtered_}
    filename_noext_multicov=${filename_compl_multicov2%.*}
    

    for complete_path_bam in $dedu_file; do
        filename_compl_bam=${complete_path_bam##*/}
        filename_noext_bam=${filename_compl_bam%.*}
        if [ "$filename_noext_bam" == "$filename_noext_multicov" ] ; then
            echo "normalizing $filenames_compl_multicov"
            echo "SAMTOOLS ---> $filename_compl_bam"

            #attivo ambiente conda samtools
            source /home/stefano/miniconda3/etc/profile.d/conda.sh
            conda activate samtools

            Rscript $script_path $complete_path_multicovs $out_bed $complete_path_bam
            echo "done" ; echo ; echo
        fi
    done
done


 
echo;echo;echo;echo;echo;echo

#--------------------------------------------------------------------------------------------------------------------------------
#|                                                                                                                              |
#|                                                         RANKING INDEX                                                        |
#|                                                                                                                              |
#--------------------------------------------------------------------------------------------------------------------------------
mkdir "$2/RI_peak"



source /home/stefano/miniconda3/etc/profile.d/conda.sh
conda deactivate


echo "-------------------------------------------------------------------------------------------------------------------------------"
echo "|                                                                                                                             |"
echo "|                                                      RANKING INDEX                                                          |"
echo "|                                                                                                                             |"
echo "-------------------------------------------------------------------------------------------------------------------------------"


#input ranking path
norm_file="$2/nscore_peak/*.bed"

#ranking index output
out="$2/RI_peak"

script_path="/home/stefano/Documents/IFO/folgiero/script/sirianalysis/RI.R"

for complete_path in $norm_file; do
    echo "RANKING $filename_1"
    Rscript $script_path $complete_path $out
    echo "done!" ; echo
done


echo;echo;echo;echo;echo;echo

#--------------------------------------------------------------------------------------------------------------------------------
#|                                                                                                                              |
#|                                                         SHARING INDEX                                                        |
#|                                                                                                                              |
#--------------------------------------------------------------------------------------------------------------------------------


#if option -v is present skip all the sharing index calc and merging
if [ $SHARING == 1 ]; then

echo "-------------------------------------------------------------------------------------------------------------------------------"
echo "|                                                                                                                             |"
echo "|                                                      SHARING INDEX                                                          |"
echo "|                                                                                                                             |"
echo "-------------------------------------------------------------------------------------------------------------------------------"


#intersect all peaks
mkdir "$2/multiinter"

bedfiles="$2/filtered_peak/*.bed"
outmultiinter="$2/multiinter"

source /home/stefano/miniconda3/etc/profile.d/conda.sh
conda activate ifo-base

bedtools multiinter -header -i $bedfiles > $outmultiinter/table_multiintersected.bed
echo "done!"

#Remove peaks duplicates
cat $outmultiinter/table_multiintersected.bed | uniq > $outmultiinter/unique_table_multiintersected.bed

#Full dataset merged
source /home/stefano/miniconda3/etc/profile.d/conda.sh
conda activate ifo-base


bedtools merge -i $outmultiinter/unique_table_multiintersected.bed -c 5 -o distinct > $outmultiinter/SI_for_peaks.bed

rm "$outmultiinter/table_multiintersected.bed" "$outmultiinter/unique_table_multiintersected.bed"


#call python script that calculate number of samples which condivide same peak
python3 $python_path $outmultiinter/SI_for_peaks.bed $outmultiinter/SI_for_peaks.bed

#--------------------------------------------------------------------------------------------------------------------------------
#|                                                                                                                              |
#|                                                          MERGE SI RI                                                         |
#|                                                                                                                              |
#--------------------------------------------------------------------------------------------------------------------------------

mkdir "$2/SI_RI"
cp "$2/multiinter/SI_for_peaks.bed" "$2/SI_RI/SI_RI_tab.bed"

echo "-------------------------------------------------------------------------------------------------------------------------------"
echo "|                                                                                                                             |"
echo "|                                                       MERGE SI RI                                                           |"
echo "|                                                                                                                             |"
echo "-------------------------------------------------------------------------------------------------------------------------------"


SI_RI_tab="$2/SI_RI/SI_RI_tab.bed"
SI_RI_tab_sed="$2/SI_RI/SI_RI_tab_sed.bed"
RI_tabs="$2/RI_peak/*.bed"
out="$2/SI_RI"

HEAD="chro"$'\t'"start"$'\t'"end"$'\t'"SI"

source /home/stefano/miniconda3/etc/profile.d/conda.sh
conda activate ifo-base


for complete_path in $RI_tabs; do 
    filename=${complete_path##*/}
    filename_noext=${filename%.*.*}
    filename_temp="temp_$filename"
    
    HEAD+=$'\t'"${filename_noext}"
    
    sort -u -k1,1 -k2,2n $complete_path | cut -d$'\t' -f 1,2,3,14 --output-delimiter=$'\t' > $out/$filename_temp
    
    #optimum = f 0.4
    bedtools intersect -a $SI_RI_tab -b $out/$filename_temp -wa -wb -loj -f 0.4 | awk '{$(NF-1)=$(NF-2)=$(NF-3)=""};1' | uniq | sed 's/    /\t/g' | sed 's/ /\t/g' | sponge $SI_RI_tab
    rm $out/$filename_temp
    echo "merged $filenames"
done
sed -i '1 i\'"$HEAD" $SI_RI_tab
sed 's/\./0/g' $SI_RI_tab >  $SI_RI_tab_sed

fi
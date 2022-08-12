home_dir="/home/fridald4/projects/def-gsarah/fridald4/hematuria_project"

# create empty table with header to put in all results (TISSUE, SYMBOL, ZSCORE, PVALUE)
awk 'NR==1' ${home_dir}/spredixcan_output/eqtl/hematuria_TOPMed_imputed_eqtl_Kidney_Cortex.csv | awk '{print "tissue",$2,$3,$5}' > ${home_dir}/comparison_multixcan_predixcan_results_table.txt

for tissue in Kidney_Cortex Pancreas Esophagus_Mucosa Skin_Not_Sun_Exposed_Suprapubic Skin_Sun_Exposed_Lower_leg
do
	grep -Ew 'PDPN|CCDC97' ${home_dir}/spredixcan_output/eqtl/hematuria_TOPMed_imputed_eqtl_${tissue}.csv | awk -v var1=${tissue} '{print var1,$2,$3,$5}' >> ${home_dir}/comparison_multixcan_predixcan_results_table.txt 
 	grep -Ew 'intron_1_13584341_13607173|intron_2_227078077_227080443|intron_12_69570517_69571276|intron_16_58040698_58041617|intron_19_41320470_41322595|intron_22_50626278_50626591' ${home_dir}/spredixcan_output/sqtl/hematuria_TOPMed_imputed_sqtl_${tissue}.csv | awk -v var1=${tissue} '{print var1,$2,$3,$5}' >> ${home_dir}/comparison_multixcan_predixcan_results_table.txt
done

grep -Ew 'PDPN|CCDC97' ${home_dir}/multixcan_output/hematuria_eqtl_smultixcan.txt | awk '{print "multixcan",$2,$15,$3}' >> ${home_dir}/comparison_multixcan_predixcan_results_table.txt
grep -Ew 'intron_1_13584341_13607173|intron_2_227078077_227080443|intron_12_69570517_69571276|intron_16_58040698_58041617|intron_19_41320470_41322595|intron_22_50626278_50626591' ${home_dir}/multixcan_output/hematuria_sqtl_smultixcan.txt | awk '{print "multixcan",$2,$15,$3}' >> ${home_dir}/comparison_multixcan_predixcan_results_table.txt


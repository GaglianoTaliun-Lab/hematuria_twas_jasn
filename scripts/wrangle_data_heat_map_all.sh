#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=1:00:00
#SBATCH --array=1-49
#SBATCH --job-name=wrangle_for_heat_map
#SBATCH --output=slurm-%x-%a.out
#SBATCH --error=slurm-%x-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G

HOMEPATH="/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project"
TISSUE=$(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/fridald4/hematuria_scripts/list_of_tissues.txt)

# create new files to print only tissue, chr, pos, gene name, zscore, pvalues, effect size and SE from annotated outputs
# these outputs can be afterwards merged by rows and used for heatmap and forest plots
awk -v tissue="$TISSUE" '{OFS = "\t"} NR==1{print "tissue", $2, "chromosome", "position", $3, $5, "effect_size"} NR>1{print tissue, $2, $11, $12, $3, $5, $4}' $HOMEPATH/spredixcan_output/gene_mapped_results/mapped_hematuria_TOPMed_imputed_eqtl_${TISSUE}.csv > $HOMEPATH/spredixcan_output/heatmap_data/hematuria_TOPMed_imputed_eqtl_all_${TISSUE}.out
awk -v tissue="$TISSUE" '{OFS = "\t"} NR==1{print "tissue", "gene_name", "chromosome", "position", "zscore", "pvalue", "effect_size"} NR>1{print tissue, $1, $2, $3, $5, $7, $8}' $HOMEPATH/spredixcan_output/gene_mapped_results/hematuria_TOPMed_imputed_sqtl_${TISSUE}.tsv > $HOMEPATH/spredixcan_output/heatmap_data/hematuria_TOPMed_imputed_sqtl_all_${TISSUE}.out

# afterwards, merge all the tissues into one file:
# HOMEPATH="/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project"
# awk 'NR==1' $HOMEPATH/spredixcan_output/heatmap_data/hematuria_TOPMed_imputed_eqtl_all_Whole_Blood.out > $HOMEPATH/spredixcan_output/heatmap_data/all_genes_eqtl_for_foresplot.out
# awk 'NR==1' $HOMEPATH/spredixcan_output/heatmap_data/hematuria_TOPMed_imputed_sqtl_all_Whole_Blood.out > $HOMEPATH/spredixcan_output/heatmap_data/all_genes_sqtl_for_foresplot.out
# for i in {1..49} ; do TISSUE=$(sed -n ${i}p /scratch/fridald4/hematuria_scripts/list_of_tissues.txt) ; awk 'NR>1' $HOMEPATH/spredixcan_output/heatmap_data/hematuria_TOPMed_imputed_eqtl_all_${TISSUE}.out >> $HOMEPATH/spredixcan_output/heatmap_data/all_genes_eqtl_for_foresplot.out ; awk 'NR>1' $HOMEPATH/spredixcan_output/heatmap_data/hematuria_TOPMed_imputed_sqtl_all_${TISSUE}.out >> $HOMEPATH/spredixcan_output/heatmap_data/all_genes_sqtl_for_foresplot.out ; done

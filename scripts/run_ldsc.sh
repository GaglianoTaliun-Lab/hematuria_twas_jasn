#!/bin/bash
#SBATCH --account=def-gsarah
#SBATCH --time=1:00:00
#SBATCH --array=1-7
#SBATCH --job-name=ldsc
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G

module load StdEnv/2020
module load python/2.7
module load scipy-stack/2020a

virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --upgrade pip
# pip install bitarray
# pip install pybedtools==0.7.10
pip install contextlib2==0.2
pip install importlib_metadata==1.5.0
pip install jsonschema==2.6.0
pip install six==1.16.0+computecanada
pip install bitarray==0.8.0

project_dir="/home/fridald4/projects/def-gsarah/fridald4/renal_genetics_project"

pip install -r ${project_dir}/ldsc/CC_requirements.txt

pheno1=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/MR_analysis/ldsc_corr/list_traits.txt)
n_studies=$(wc -l ${project_dir}/MR_analysis/ldsc_corr/list_traits.txt | awk '{print $1}')

for ((i=1;i<=n_studies;i++))
do
	pheno2=$(sed -n ${i}p ${project_dir}/MR_analysis/ldsc_corr/list_traits.txt)
	python ${project_dir}/ldsc/ldsc.py \
		--ref-ld-chr ${project_dir}/reference_data/eur_w_ld_chr/ \
		--out ${project_dir}/ldsc_corr/rg_output/rg_${pheno1}_${pheno2} \
		--rg ${project_dir}/MR_analysis/ldsc_corr/${pheno1}.sumstats.gz,${project_dir}/MR_analysis/ldsc_corr/${pheno2}.sumstats.gz \
		--w-ld-chr ${project_dir}/reference_data/eur_w_ld_chr/ 
done

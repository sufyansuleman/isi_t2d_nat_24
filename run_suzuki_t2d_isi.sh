#!/bin/bash
#SBATCH --job-name=t2d_assoc
#SBATCH --output=t2d.out
#SBATCH --error=t2d.err
#SBATCH --time=10:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=1

# 🔥 RESET ENVIRONMENT (CRITICAL)
module purge

# initialize modules if needed (safe to include)
source /usr/share/Modules/init/bash

# load compatible stack
module load gcc/11.2.0
module load openjdk/20.0.0
module load R/4.3.3

# debug (do NOT remove for now)
which gcc
which R
which Rscript
R --version

cd /projects/glostrup-AUDIT/people/rnh585/insulin_sensitivity/gwas/scripts/

Rscript suzuki_t2d_isi_assoc.R


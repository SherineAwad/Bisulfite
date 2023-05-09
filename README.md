
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)


# Bisulfite Sequencing pipeline

Then run: snakemake -jnumber_of_cores, for example for 5 cores use:

     snakemake -j5 


and for a dry run use:

    snakemake -j1 -n 



and to print the commands in a dry run use:

      snakemake -j1 -n -p 



To use another config file use:

     snakemake -j1 -p --configfile configfilehere.yaml


For the sake reproducibility, use conda to pull same versions of tools. Snakemake and conda have to be installed in your system:

     snakemake --cores --use-conda

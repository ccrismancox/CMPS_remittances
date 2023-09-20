# Replication package for: "Remittances, terrorism, and democracy"

## Overview
This replication archive contains the data and code for replicating the article  "Remittances, terrorism, and democracy." The code was run using R 4.3.1 ("Beagle Scouts") on a Ubuntu 22.04 computer with the following session information

    R version 4.3.1 (2023-06-16)
    Platform: x86_64-pc-linux-gnu (64-bit)
    Running under: Ubuntu 22.04.3 LTS
    
    Matrix products: default
    BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0
    
    locale:
    [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
    [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
    [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    
    time zone: America/Chicago
    tzcode source: system (glibc)
    
    
    attached base packages:
    [1] splines   stats     graphics  grDevices utils     datasets  methods   base     
    
    other attached packages:
     [1] ggplot2_3.4.3      matrixStats_1.0.0  margins_0.3.26     lme4_1.1-34        Matrix_1.6-0      
     [6] pscl_1.5.5.1       stargazer_5.2.3    stringr_1.5.0      xtable_1.8-4       lmtest_0.9-40     
    [11] zoo_1.8-12         maxLik_1.5-2       miscTools_0.6-26   nonnest2_0.5-6     modelsummary_1.4.2
    [16] MASS_7.3-60        sandwich_3.0-2     car_3.1-2          carData_3.0-5      data.table_1.14.8 




## Main level

The main level contains three files and three folders. The files are:

-   `README.md` This file the in text format
-   `master.R` An R code file that replicates the paper. This file changes the working directory to the `code` folder and runs all of the R scripts. Outputs to `output/masterRun.txt`


The three folders are described below

## Data

The folder `data` contains 2 files:

-   `mainRemittanceData.rdata` This is the main data set used in the analysis.
-   `oecdmembers.csv` List OECD members

### Code

The folder `code` contains 12 R scripts

-   `packages.r` Installs packages and details the version information as of 9/20/2023
-   `nbreg.r` This file contains helper functions for fitting the negative binomials when glm.nb issues warnings
-   `stargazerNoteCorrection.r` code to fix a bug in stargazer's note options
-   `summary.R` Produces Tables 1, 2, and A1
-   `main_analysis.R` Fits the main models and produces `output/mainResultsRemittances.rdata`
-   `main_table_figures.R` Uses the results in `output/mainResultsRemittances.rdata` to produce Tables 3, 4, & B1 and Figures 1 and C1
-   `vdem_partof_table5.R` Fits models 7-8 and saves them in `output/vdem_competition.rdata`
-   `table5.R` Fits models 5-6 and uses `output/vdem_competition.rdata` to produce Table 5
-   `table6.R` Produces Table 6
-   `appendixD1.R` Produces Table D1
-   `appendixD2.R` Produces Table D2
-   `appendixD3.R` Produces Tables D3 and D4
-   `appendixD5.R` Produces Table D5

## Output

The folder `output` contains 5 files. The generation and contents of these files are described above

## Running the code

Files may be run individually from the `code` folder or the file `master.r` can be run from the main folder.

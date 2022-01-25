# Sudo-Tanaka-2021-JMCB_ReplicationCodes
This is the replication codes for Sudo and Tanaka (2021), "Quantifying Stock and Flow Effects of QE," Journal of Money Credit and Banking. The codes are identical to those available on the journal website.


## A) Requirements
 - MATLAB
 - DYNARE
 - 00_tool (contained in this zip file)

  * Please add `00_tools` to your MATLAB search path, in advance of running the replication codes.

  * We used MATLAB 2016a and DYNARE 4.4.3.

  * Older versions of DYNARE (including v4.4.3) are available on the official webpage of DYNARE ( https://www.dynare.org/release/ ).



## B) Description of Replication Codes

### B-0) Overview 

[../00_tools]
 This directory contains MATLAB functions which are required to compute and 
 plot impulse response functions.

[../01_Baseline]
 This directory contains MATLAB and DYNARE codes for our "baseline" model.
 The codes in this directory are required to replicate the following
 tables and figures.
     Table 3, Table 4, Figure 4, Figure 5, Figure 6(1)

[../02_Alternative1]
 This directory contains MATLAB and DYNARE codes for our "Alternative 1" 
 model. The codes in this directory are required to replicate the following
 table and figures.
     Table 4, Figure 6(2), Figure 7

[../03_Alternative2]
 This directory contains MATLAB and DYNARE codes for our "Alternative 2" 
 model. The codes in this directory are required to replicate the following
 table and figures.
     Table 4, Figure 6(3), Figure 7

[../04_Alternative3]
 This directory contains MATLAB and DYNARE codes for our "Alternative 3" 
 model. The codes in this directory are required to replicate the following
 table and figures.
     Table 4, Figure 6(4), Figure 7


### B-1) Replicate Table 3

1. Run `Step1_estimation.m` in 01_Baseline.
    
2. Open `Baseline_estimation_results.mat,` which DYNARE automatically will make in the estimation procedure.
    
    oo_.posterior_mean.parameters:
        Posterior means of the estimated parameters

    oo_.posterior_hpdsup.parameters
        Upper bounds of the posterior 90% intervals

    oo_.posterior_hpdinf.parameters
        Lower bounds of the posterior 90% intervals

### B-2) Replicate Figure 4

[IRFs with commitment to ZLB, depicted as thick lines in Figure 4]
 1. Run `Step2_IRF_with_commitment.m` in 01_Baseline.

 2. IRFs will appear in a new graph window.


[IRFs without commitment to ZLB, depicted as thin lines in Figure 4]
 1. Run `Step3_IRF_without_commitment.m` in 01_Baseline.

 2. IRFs will appear in a new graph window.


### B-3) Replicate Figure 5

Note: This step takes much time (approximately one hour), because DYNARE is iteratively implemented using 1000 sets of drawn parameters.

 1. Run `Step4_historical_decomp.m` in 01_Baseline.

 2. The results will be stored in `Historical_Decomposition_Baseline.xlsx.`


### B-4) Table 4

1. For each alternative model, run `Step1_estimation.m` in the corresponding directory.

2. Open `AlternativeX_estimation_results.mat` (X = 1, 2, 3)
    
    oo_.posterior_mean.parameters:
        Posterior means of the estimated parameters

    oo_.posterior_hpdsup.parameters:
        Upper bounds of the posterior 90% intervals

    oo_.posterior_hpdinf.parameters:
        Lower bounds of the posterior 90% intervals


### B-5) Figure 6

 1. For each alternative model, run `Step2_IRF_without_commitment.m` in the corresponding directory.

 2. IRFs will appear in a new graph window.


### B-6) Figure 7

Note: This step takes much time (approximately one hour), because DYNARE is iteratively implemented using 1000 sets of drawn parameters.

 1. For each alternative model, run `Step3_historical_decomp.m` in the corresponding directory.

 2. The results will be stored in `Historical_Decomposition_AltX.xlsx.` (X = 1, 2, 3)



## C) References

Adjemian, S., H. Bastani, M. Juillard, F. Karame, J. Maih, F. Mihoubi, W. Mutschler, G. Perendia, J. Pfeifer, M. Ratto and S. Villemot, "Dynare: Reference Manual, Version 4," 2011, CEPREMAP.

   * See the manual for the software requirements, installation, and usage of DYNARE.


Chen, H., V. Curdia, and A. Ferrero, "Technical Appendix to the Macroeconomic Effects of Large-Scale Asset Purchase Programmes," 2012, Wiley.
   https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fj.1468-0297.2012.02549.x&file=ecoj2549_sm_Appendix.pdf
     
   * Our simulation method used in Step B-2 and B-5 are based on their Appendix F.2.


---
Copyright 2018-2020 by Nao Sudo and Masaki Tanaka (Bank of Japan)

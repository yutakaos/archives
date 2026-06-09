# Information for the compiled data
Ushio (2022) "Interaction capacity as a potential driver of community diversity"


## "data_sample.csv" includes following information.
- Sample_ID: Sample identifier
- Sample_Name2: Sample name for the analytical purpose (identical to Sample ID)
- Description: More informative sample identifier.
- date: Sampling date
- plot: Sampling plot
- sample_nc: "sample" or "nc" (negative control, not included in the compiled data)
- filt045_ml: Filtered water volume when using 0.45-um filter cartridge
- filt022_ml: Filtered water volume when using 0.22-um filter cartridge

In the following columns, "XXXX" indicates "Prok" (16S), "Fungi" (ITS), "Inv" (COI), or "Euk" (18S).
- XXXX_dada2_input: Sequence reads of XXXX region that were inputted to DADA2
- XXXX_dada2_prop: (DADA2 processed sequence reads)/(Input sequence reads)
- XXXX_STD_all: Sequence reads of standard DNA
- XXXX_NonSTD_all: Sequence reads of non-standard DNA
- XXXX_STD_prop: Proportion of standard DNA sequence reads
- XXXX_STD_r2: R2 of the linear regression line of standard DNA copy numbers v.s. sequence reads
- XXXX_STD_coef: Coefficient of the linear regression line of standard DNA copy numbers v.s. sequence reads 
- XXXX_dna_copy_sum: Total DNA copy number of XXXX region


## "data_asvtable.csv" includes DNA copy numbers / ml water in each sample.
- Columns are taxa IDs and rows are sample IDs.


## "data_taxa.csv" inculudes following information.
- Taxa_ID: Taxa indentifier
- seq: Sequence of the ASV
- seqlen: Sequence length of the ASV
- entropy: Entropy of the time series data of the ASV (this was used to remove non-informative ASVs)
- miseq_run: MiSeq run ID to generate the ASV


## "data_climate.csv" inculudes following information.
- date: Date
- plot: Plot
- temp_mean: Mean daily air temperature
- relhm_mean: Mean relative humidity
- satdf_mean: Mean saturation defficiency
- actvp_mean: Mean actual vapor pressure
- light_mean: Mean light intensity
- temp_max: Maximum daily air temperature
- relhm_max: Maximum relative humidity
- light_max: Maximum light intensity
- temp_min: Minimum daily air temperature
- relhm_min: Minimum relative humidity
- light_min: Minimum light intensity



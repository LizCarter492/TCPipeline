## Code and data for:
### Smith & Carter, **Increasing hazardous material pipeline failures during tropical cyclones: trends and future implications**, ERL, *in submission*
Tropical cyclones (TCs, including hurricanes) have caused nearly $1.5 trillion in damages in the US since 1980—more than all other kinds of natural disasters combined. TCs frequently intersect with America’s hazardous materials pipeline network: a vast network of pipes that transport the majority of crude oil and natural gas to the public. Since 1970, over 30,000 pipeline failures have been reported, causing hundreds of deaths and billions in damages. TCs are under-recognized as causes of pipeline failures, mainly because of how failures are reported. In this study, we make a conservative inventory of pipeline failures that are associated with TCs. TCs impacting pipeline infrastructure have intensified since 1970 due to climate change. Because of this, the frequency of storm-associated pipeline failures has doubled since 1975. If current trends continue, the number of pipeline failures per TC could increase by 50-600%, putting more pressure on the already vulnerable pipeline infrastructure. Spatial variability in the frequency of pipeline failures after TCs gives some clues as to causes, including coastal storm surges, tornados, and mechanical forcing from shrink swell clays. These findings highlight the urgent need for TC-specific regulatory measures to protect hazardous material pipelines from intensifying TCs.

<p align="center">
<img src="https://github.com/LizCarter492/TCPipeline/blob/main/figures/Figure1.png" width="500" >
</p>
<em>Figure 1: Six regions where HMP infrastructure (major crude oil and natural gas transmission pipelines, EIA 2020) intersect with tropical cyclone storm tracks in the continental United States. For each region, table indicates total 1970-2022 TC-associated failures, maximum TC windspeed, length of crude oil pipeline, and length of natural gas pipeline</em>.


#### Use this repository to replicate our analysis:

1. code:
    * 1_formatData.ipynb : reads from **raw_data** (HURDAT2 and PHMSA failure datasets), merges on spatiotemporal concomittence, saves to **analysis_data** 
    * 2_modelTrain.R : reads from **analysis_data**, fits pca, binomial, and poisson regression, saves to **model_output**
    * 3_figures.R: reads from  **raw_data, analysis_data, model_output**; generates figures, saves to **figures**  
3. raw_data:
    * NOAAHURDAT.csv : dowloaded from [insert link here]()
    * PHMSA.csv : downloaded from [insert link here]()
5. analysis_data:
    * pipeline_failures_storm.csv : full PHMSA failure dataset with information about storm overlaps
    * analysis_ready_binomial.csv : PHMSA dataset truncated to TC storm track regions, all HMP failure/TC details, with pressence/absense of failure included in "failTF" column.
    * analysis_read.csv: analysis_ready_binomial flattened on "StormID", with failTF_sum indicating the storm-level failure frequency (fs) 
7. model_output:
    * contains all data and model objects generated by 2_modelTrain.R, as .RData files
9. figures:
    * contains all figures generated by 3_figures.R

Questions/comments/pull requests? Contact ekcarter@syr.edu

   [This is a GRRIEn analysis](https://journals.ametsoc.org/view/journals/aies/2/2/AIES-D-22-0065.1.xml) 


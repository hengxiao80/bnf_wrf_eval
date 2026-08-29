As a first step, I want to use /ccsopen/home/hengxiao80/bnf_wrf_eval/satoshi_testruns/20250917lassobnfwrfhrrr3/rund1 as an example simulation to construct and test scripts needed to plot comparision plots between the simulation and the HRRR data.

Later on, I want to expand to other simulations and other analyses (ERA5) and obs. datasets.

The first variable I want to compare is radar reflectivity and OLR.

As the first step of the first step of the first step, can you check whether OLR and some kind of radar reflectivity variable are avaliable in the WRF output as well as in the HRRR data we downloaded?

Use the HRRR native levels data at satoshi_forcing_data/hrrr/hrrrnat_data. And also use analyses for now. We may want to switch to using forecast data later.


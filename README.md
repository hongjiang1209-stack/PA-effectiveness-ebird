# PA-effectiveness-ebird

This folder contains the core R scripts used in the analyses for:

“Only half of protected areas effectively conserve birds in the Contiguous United States”

Scripts are numbered approximately according to the analytical workflow.

Analyses were conducted using:

-R 4.3.2
-ArcGIS 10.8
-ArcGIS Pro

Some scripts were run on the University of York Viking2 HPC system.

The analyses used publicly available datasets, including:

-eBird Basic Dataset (EBD)
-World Database on Protected Areas (WDPA)
-MODIS net primary production data
-Global waterbody datasets
-GEBCO elevation data

File paths may require modification before execution.


Scripts included in this folder:

1.read_EBD.R
2.Merge_filter_checklists.R
3.Number_of_checklist_in_each_PA.R
4.Read_variables.R
5.Subset_checklists.R
6.Matching_per_PA.R
7.SR_estimation.R
8.Effectiveness.R
9.Complementarity.R

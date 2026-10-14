# "Integrating cell-state and enhancer constraints discriminates models of PU.1 regulation during hematopoietic differentiation"
Model and figure generating scripts for the pre-print publication. The README is structured as follows:
### folderName
1. whichFigureTheScriptIsFor: MATLAB/FILE/TO/RUN
2. anotherFigure: ANOTHER/MATLAB/FILE/TO/RUN

---
## Figure1_(quant_comparison)
1. Figure 1B and 1C: Figure1_(quant_comparison)/Fig4_replicate_dimensionless.m

## Figure1_(cell_cycle_contradiction)
2. Figure 1D: Figure1_(quant_comparison)/Plot_cellcycle_length_KuehModel_vs_Real_v3.m

## Figure23_(phen_results)
1. Figure 2BCD and Figure S1: Figure2(PU1_alone_no data_use_ps)/PU1_alone_Hill_use_ps.m
2. Figure 3BC: Figure3(PU1+CE+nofeedback_X=const)/PU1_CEBPa_Hill_noData_ps_to_find_sol_ODE_Xconst.m
3. Figure 3E: Figure3(PU1+CE+feedback+CEdilution)/PU1_CEBPa_Hill_noData_wFeedback_ps_to_find_sol_v2.m

## Figure45_(biophy_results)
1. Figure 4CD: cebpa_wPCfeedback_noAUTOfeedback_extraPenal(model1)/ODE_Local_Post_Analysis.m
2. Figure 4EF, Figure S3BCD, and Figure S4BC: debug_why_cebpa_fails/analyze_nofeedback_results.m
3. Figure S5BC: cerun_wPCfeedback_noAUTOfeedback_extraPenal(model3)/ODE_Local_Post_Analysis.m
4. Figure 5BCD: BifurcationAnalysisofModel3_noMATCONT_files/Best_pars_obtained_from_bifur.m
5. Figure 5E: NullclineForEachModel/Model1/model1_prolow_reversexy.m
6. Figure 5F: NullclineForEachModel/Model3/model3_prolow_reversexy.m
7. Figure S6BCD: cerun_wPCfeedback_wAUTOfeedback_extraPenal(model4)/ODE_Local_Post_Analysis.m
8. Figure S7BC: debug_why_cerun_noCEAUTO_fails_v2/predict_debug_why_CEAUTO_necessary_v2.m

*Note: best parameters of the "C/EBPa model" in main text are found through bifurcation analysis. Due to the enormous size of code packages needed for MATCONT, we excluded these files during upload. For any interest in re-running the bifurcation analyses in folder "BifurcationAnalysisofModel3_noMATCONT_files", simply copy and paste MATCONT packages into the folder.* 


## Error_propagation_and_normalization
1. two-dimensional uncertainty in joint measurements of PU.1 abundance and synthesis rate after being normalized to macrophage: Error_propagation_and_normalization/error_prop_concAndProdRate.m
2. uncertainty in PU.1 dilution rate converted from cell cycle length after being normalized to macrophage: Error_propagation_and_normalization/error_prop_growthRate.m
3. uncertainty in C/EBPa dilution rate converted from cell cycle length after being normalized to macrophage: Error_propagation_and_normalization/error_prop_RemRate_CEBPA.m

*Note: for reproducibility, all results in the main text are developed using the results saved in "error_prop_0216.mat"*
	  

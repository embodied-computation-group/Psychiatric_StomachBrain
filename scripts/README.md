# Scripts for psychiatric stomach-brain coupling CCA

1. Prepare stomach-brain coupling (phase locking value) data for CCA (empirical - chance PLV, and DiFuMo parcellate):
    Psychiatric_StomachBrain/scripts/CCA/CCAprepStomachBrain.ipynb

2. Create matched input matrices for CCA (X = stomach-brain coupling, Y = psych data, C = confounds. Outliers removes & rows across matrices matched as same participants):
    Psychiatric_StomachBrain/scripts/CCA/CCAprepAllData.m

3. X.mat, Y.mat & C.mat from (2.) saved in: 'Psychiatric_StomachBrain/scripts/CCA/cca_pls_toolkit-master/_Project_StomachBrain/data', and framework folder created: '/home/leah/Git/Psychiatric_StomachBrain/scripts/CCA/cca_pls_toolkit-master/_Project_StomachBrain/framework'

4. Run CCA with CCA/PLS toolkit (run on cluster via Psychiatric_StomachBrain/scripts/CCA/cca_pls_toolkit-master/cca_jobs_slurm.sh):
    Psychiatric_StomachBrain/scripts/CCA/cca_pls_toolkit-master/RunCCA.m

5. Plot CCA result (variate scatterplot, psych loading barplot, extraction of cca mode result details):
    Psychiatric_StomachBrain/scripts/CCA/cca_pls_toolkit-master/CCA_plots.m

6. Plot CCA result: stomach-brain coupling loadings projected on DiFuMo parcellated brain:
    Psychiatric_StomachBrain/scripts/CCA/CCA_plotting/CCA_Brain_Figures.ipynb

7. Create averaged summary plots of CCA loadings (for both psych and stomach-brain):
    Psychiatric_StomachBrain/scripts/CCA/CCA_plotting/Averaged_loadings_plots.ipynb

--------------------------------------------------------------------------------------------------------
Extra scripts for control analyses:

8. Prepare control standard deviation of BOLD activity for CCA (DiFuMo parcellate):
    Psychiatric_StomachBrain/scripts/CCA/Control_analyses/ControlCCAs/CCAprepControlSTD.ipynb

9. Create matched input matrices for CCA (X = Control STD BOLD, Y = psych data, C = confounds. Rows across matrices matched as same participants):
    Psychiatric_StomachBrain/scripts/CCA/Control_analyses/ControlCCAs/CCAprepAllData_controls.m

10. Repeat 3-4 but with results from (9.).

11. Run 9-10 for resting connectivity control CCA (commented in 9. replacing STD BOLD for resting connectvity data).

12. 
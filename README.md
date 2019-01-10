# GLM-CFC
A GLM based approach to assess CFC.

Code to measure cross-frequency coupling between two signals as described in this manuscript: 
*A statistical modeling framework to assess cross-frequency coupling while accounting for confounding effects*, J. Nadalin, 
L-E Martinet, A. S. Widge, S. S. Cash, U. T. Eden, M. A. Kramer, 2018.

`ExampleCode.m`: Run the cells in this file to produce example voltage traces and surfaces in (Phi_low, A_low, A_high)-space (as in Figure 4 of the manuscript). Four simulations are present: (i) no CFC, (ii) PAC only, (iii) AAC only, and (iv) both PAC and AAC.

`simfun.m` : code to simulate signals V_low, V_high, with induced cross-frequency coupling and measure output statistics R_PAC, R_AAC, and R_CFC, along with confidence intervals and p-values

`glmfun.m`: code to evaluate the coupling statistics R_PAC, R_AAC, and R_CFC, along with confidence intervals and p-values, between two signals.

`ExampleCode`: run to get example voltage traces and surfaces in Phi_low, A_low, A_high space (as in Figure 4). Four simulations are present: one with no CFC, one with PAC, one with AAC, and one with both PAC and AAC

`glmfun_with_indicator`: an example of how to update glmfun to test for effect of condition (e.g. pre and post stimuli) on coupling

Any questions/comments please direct to Jessica Nadalin (jnadalin@bu.edu) or Mark Kramer (mak@bu.edu)

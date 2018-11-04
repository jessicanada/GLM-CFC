# GLM-CFC
A GLM based approach to measure CFC

Code to measure cross-frequency coupling between two signals as described in this manuscript: [A statistical modeling framework to assess cross-frequency coupling while accounting for confounding effects]

simfun: code to simulate signals V_low, V_high, with induced cross-frequency coupling and measure output statistics R_PAC, R_AAC, and R_CFC, along with confidence intervals and p-values

glmfun: code to evaluate coupling statistics R_PAC, R_AAC, and R_CFC, along with confidence intervals and p-values, between two signals

ExampleCode: run to get example voltage traces and surfaces in Phi_low, A_low, A_high space (as in Figure 4). Four simulations are present: one with no CFC, one with PAC, one with AAC, and one with both PAC and AAC

Any questions/comments please direct to Jessica Nadalin (jnadalin@bu.edu)

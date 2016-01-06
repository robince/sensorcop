
This repository contains the code for the analysis and simulations in:

RAA Ince, BL Giordano, C Kayser, GA Rousselet, J Grouss and PG Schyns  
"A statistical framework based on a novel mutual information estimator utilizing a Gaussian copula"

This is the full history of the active research code. It is unlikely to run unaltered on a fresh machine without changing data and toolbox paths as required. There are also several mex functions which need compilation.

The `eegface` folder contains the code used for the discrete stimulus event-related EEG example in Section 4.1. 

The `megspeech` folder contains the code used for the continuous stimulus continuous MEG example in Section 4.2.

The particular scripts used for each figure are:

- **Figure 1** : [`fig_entropy_examples.m`](fig_entropy_examples.m)
- **Figure 2** : `info_scatters.py`
- **Figure 3** : `corr_vs_info.m`
- **Figure 4A** : `eegface/fig_event_related_design.m`
- **Figure 4(B,D)** : `megspeech/fig_conf_design.m`
- **Figure 5** : `fig_cmi_example.m`
- **Figure (8,9)** : `fig_copula_examples.m`


This repository contains the code for the analysis and simulations in:

RAA Ince, BL Giordano, C Kayser, GA Rousselet, J Gross and PG Schyns  
"A statistical framework for neuroimaging data analysis based on mutual information estimated via a Gaussian copula"  
[bioRxiv:](http://biorxiv.org/content/early/2016/03/16/043745) [doi:10.1101/043745](http://dx.doi.org/10.1101/043745)

This is the full history of the active research code. It is unlikely to run unaltered on a fresh machine without changing data and toolbox paths as required. There are also several mex functions which need compilation (see [`make.m`](make.m)).

The raw data and plotted figure data will be available in the Dryad data package:

RAA Ince, BL Giordano, C Kayser, GA Rousselet, J Gross and PG Schyns  
"Data from: A statistical framework based on a novel mutual information estimator utilizing a Gaussian copula"

The public version of the Gaussian-Copula Mutual Information (GCMI) estimator is available [here](https://github.com/robince/gcmi).

The `eegface` folder contains the code used for the discrete stimulus event-related EEG example in Section 4.1. 

The `megspeech` folder contains the code used for the continuous stimulus continuous MEG example in Section 4.2.

The `eegeye` folder contains the code for the continuous stimulus event-related EEG interaction information examples in Section 4.3.

The particular scripts used for each figure are:

- **Figure 1** : [`fig_entropy_examples.m`](fig_entropy_examples.m)
- **Figure 2** : [`info_scatters.py`](info_scatters.py)
- **Figure 3** : [`corr_vs_info.m`](corr_vs_info.m)
- **Figure 4A** : [`eegface/fig_event_related_design.m`](eegface/fig_event_related_design.m)
- **Figure 4(B,D)** : [`megspeech/fig_cont_design.m`](megspeech/fig_cont_design.m)
- **Figure 5** : [`fig_cmi_example.m`](fig_cmi_example.m)
- **Figure (8,9)** : [`fig_copula_examples.m`](fig_copula_examples.m)
- **Figure 10** : [`fig_spectral_example_separate.m`](fig_spectral_example_separate.m)
- **Figure 11A** : [`eegface/plot_full_erp.m`](eegface/plot_full_erp.m)
- **Figure 11B** : [`eegface/plot_example_tc.m`](eegface/plot_example_tc.m)
- **Figure 11C** :
    + Run simulations (number of trials): [`eegface/full_erp_maxstats_1d.m`](eegface/full_erp_maxstats_1d.m)
    + Plot results (number of trials): [`eegface/plot_cm_stats.m`](eegface/plot_cm_stats.m)
    + Run simulations (robustness): [`eegface/full_erp_robustness_1d.m`](eegface/full_erp_robustness_1d.m)
    + Plot results (robustness): [`eegface/plot_robustness1d.m`](eegface/plot_robustness1d.m)
- **Figure 12A** : [`megspeech/sens_delay_different_sigs.m`](megspeech/sens_delay_different_sigs.m)
- **Figure 12B** : [`megspeech/plot_example_chan.m`](megspeech/plot_example_chan.m)
- **Figure 12C** :
    + Run simulations (data length): [`megspeech/sens_delay_maxstats_1d.m`](megspeech/sens_delay_maxstats_1d.m)
    + Plot results (data length): [`megspeech/plot_cm_stats.m`](megspeech/plot_cm_stats.m)
    + Run simulations (robustness): [`megspeech/sens_delay_robustness_1d.m`](megspeech/sens_delay_robustness_1d.m)
    + Plot results (robustness): [`megspeech/plot_robust1d.m`](megspeech/plot_robust1d.m)
- **Figure 13(A-E)** : 
    + MI calculations: [`eegeye/calc_info_temporal.m`](eegeye/calc_info_temporal.m)
    + Plot results: [`eegeye/fig_parainfo_int_gradient.m`](eegeye/fig_parainfo_int_gradient.m)
- **Figure 13F** : [`eegeye/info_accum.m`](eegeye/info_accum.m)
- **Figure 13G** : [`megspeech/plot_example_chan.m`](megspeech/plot_example_chan.m)
- **Figure 14A** : [`gausssim/gausssim_1d.m`](gausssim/gausssim_1d.m)
- **Figure 14B** : [`gausssim/gausssim_2d.m`](gausssim/gausssim_2d.m)
- **Figure 14C** : [`eegface/info_bias_2d_jack.m`](eegface/info_bias_2d_jack.m)
- **Figure 14D** : [`megspeech/meg_info_bias_2d_jack.m`](megspeech/meg_info_bias_2d_jack.m)

Questions / comments : robince@gmail.com

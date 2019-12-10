# Data and scripts for Matthewson, Woo, Reddan, & Wager (2019)

All scripts are written in matlab (version > 2016b) using matlab live script (mlx). Unfortunately, github does not support .mlx yet. To see the scripts, you need to open the scripts using matlab.

### Dependencies

- [CanlabCore](https://github.com/canlab/canlabcore)
- [CocoanCORE](https://github.com/cocoanlab/cocoanCORE)

### Reference

Gordon Matthewson\*, Choong-Wan Woo\*, Marianne C. Reddan, Tor D. Wager<sup>§</sup> (in press) Cognitive self-regulation influences pain-related physiology, _PAIN_

*co-first authors

<sup>§</sup>corresponding author: Tor D. Wager (tor.d.wager@dartmouth.edu)

PAIN: [link](https://journals.lww.com/pain/Abstract/publishahead/Cognitive_self_regulation_influences_pain_related.98669.aspx)

Preprint: [biorxiv](https://www.biorxiv.org/content/10.1101/361519v1)

Video abstract: [Youtube](https://www.youtube.com/watch?v=R1QtvyAt-F8)

[![Video abstract](https://img.youtube.com/vi/R1QtvyAt-F8/0.jpg)](https://www.youtube.com/watch?v=R1QtvyAt-F8 "Video abstract")


### Abstract

Cognitive self-regulation can shape pain experience, but its effects on affects autonomic responses to painful events is unclear. In this study, participants (N = 41) deployed a cognitive strategy based on reappraisal and imagination to regulate pain up or down on different trials while skin conductance responses (SCR) and electrocardiogram (ECG) activity were recorded. Using a machine learning approach, we first developed stimulus-locked SCR and ECG physiological markers predictive of pain ratings. The physiological markers demonstrated high sensitivity and moderate specificity in predicting pain across two datasets, including an independent test dataset (N = 84). When we tested the markers on the cognitive self-regulation data, we found that cognitive self-regulation had significant impacts on both pain ratings and pain-related physiology in accordance with regulatory goals. These findings suggest that self-regulation can impact autonomic nervous system responses to painful stimuli and provide pain-related autonomic profiles for future studies.

About the pain-related skin conductance response (SCR)/Electrodermal Response (EDR) measure:

The weight pattern to apply to new data can be found in https://github.com/canlab/cognitive_regulation_physiology/tree/master/data

SCR_weights.mat
ECG_weights.mat

Each one has the 500x1 vector for intensity and unpleasantness.

To use it, you would need to apply it to an EDA time series locked to stimulus onset. It was validated to pain intensity under variations in painful stimulation, with limited specificity to pain, as described in our paper (Matthewson, Woo et al. in press, Pain; Bioarxiv https://www.biorxiv.org/content/10.1101/361519v2.abstract). However, it has not been broadly validated to be specific to pain, so there is no set cutoff for the amplitude of the integrated measure that signals “pain”. 

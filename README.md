# ECG_project
Electrocardiograph (ECG) systems measure electrical activity of muscle contractions via carefully attached electrodes. Analog band pass filters are a standard method of signal processing to eliminate noise components from the ECG signal. This work assessed the magnitude and phase response of an analog band pass filter used in an industrial ECG design to ascertain the extent of the signal distortion. A series of all pass filters were used to decrease the phase distortion by 50 \% at a trade off of slightly more signal delay. Flatness of the group delay served as an optimization parameter to configure the all pass system.

## Approach
All pass filters can be cascaded to produce a more desirable phase response. A parameter sweep over phase angle and radius r was carried out to create a range of pole/reciprocal zero locations. Phase angle ranged from 0 to &pi; with an increment of .05 and r from 0 to 1 with an increment of .01. The geometric mean provided better results during optimization testing than the arithmetic mean. For each iteration of the parameter sweep, the geometric mean, &mu; is calculated for the cascade of the original band pass with this candidate all pass. After each sweep, the all pass filter that produced the minimum &mu; and is below a certain threshold value is added as a stage to the all pass system. The process is then repeated with this new band pass all pass transfer function. The goal is to progressively make the group delay more constant in the frequency band of interest until it meets a sufficient exit condition. The figures below shows the Bode plot for a bandpass filter cascaded with six all pass systems and the right plot shows the pole-zero plot.

<img src="img/bode_cas.png?raw=true" width="345">             <img src="img/hd_pz.jpg?raw=true" width="365">

## Results
Results of the optimization are found in the figures below. In the left figure, we see that the addition of all pass stages ceased to make the group delay more constant after six iterations. In the right figure, we see the bandpass filter cascaded with the all pass systems does not visually distort the original signal after the filtering.

<img src="img/group_delay_seq.png?raw=true" width="366">             <img src="img/hd_filt.jpg?raw=true" width="370">

The table below shows the root-mean-square error (RMSE) between the input signal and the filtered output signal. BPF denotes the standard bandpass filter and BPF-AP denotes the bandpass filter cascaded with all pass filters to linearize the phase. We see the RMSE is decreased by 50\%.

| Filter System | RMSE  |
| ------------- | ------|  
| BPF           | 0.0034|
| BPF-AP        | 0.0018|

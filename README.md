# ECG_project
Electrocardiograph (ECG) systems measure electrical activity of muscle contractions via carefully attached electrodes. Analog band pass filters are a standard method of signal processing to eliminate noise components from the ECG signal. This work assessed the magnitude and phase response of an analog band pass filter used in an industrial ECG design to ascertain the extent of the signal distortion. A series of all pass filters were used to decrease the phase distortion by \%50 at a trade off of slightly more signal delay. Flatness of the group delay served as an optimization parameter to configure the all pass system.

#Approach
All pass filters can be cascaded to produce a more desirable phase response. A parameter sweep over phase angle and radius $r$ was carried out to create a range of pole/reciprocal zero locations. Phase angle ranged from 0 to $\pi$ with an increment of .05 and $r$ from 0 to 1 with an increment of .01. The criterion for the optimal all pass locations were:

```math
    \text{min} \hspace{1mm} \left\{ \mu = \left(\prod_{k=1}^{N} |\frac{d^2}{d\theta^2}\phase{H_{BP}H_{AP}(e^{j\omega k}})|\right)^{\frac{1}{N}}  \right \}
```
 where $N$ is equivalent to 4000 points over a frequency band of 5 - 195 Hz. The geometric mean provided better results during optimization testing than the arithmetic mean. For each iteration of the parameter sweep, $\mu$ is calculated for the cascade of the original band pass with this candidate all pass. After each sweep, the all pass filter that produced the minimum $\mu$ and is below a certain threshold value is added as a stage to the all pass system. The process is then repeated with this new band pass all pass transfer function. The goal is to progressively make the group delay more constant in the frequency band of interest until it meets a sufficient exit condition. 

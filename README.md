**Here is the code for the work ``General auxiliary variables with approximate neighbourhood interference''.** R version: 4.1.2.



## Overview

We provide the experimental scripts:

- Experiments in the main body of the paper:
  - **`ex`**: Real-world experiment in our main text.
  - **`syn`**: Synthetic experiment in our main text.

- Experiments in Appendix (Additional Experimental Results):
  - **counter: `negative_result.R`**: Counter example in our Appendix S.4.
  - **modify: `linear_final_attempt_bni.R`**: bandwith selection in our Appendix S.5.
  - **modify: `syn_modify_misspecification_i.R`**: misspecification in our Appendix S.5. Note that each round, using a different random seed induces a different misspecification scenario, which in turn leads to small fluctuations in the reported table results (you can try starting it several times if any abnormalities occur during the misspecified-network-construction). Throughout this process, our method is relatively more robust to the choice of random seed.



## Running the Experiments

Please run the following commands to execute the experiments. Please make sure to change the data storage path to the one you need.


```bash
# Example:

# Run the counter example
Rscript negative_result.R

# Run the synthetic experiments
Rscript linear_final_attempt.R
Rscript nonlinear_final_attempt.R
```


## License

This code is for academic and research purposes.

<!--




## Synthetic result:


| Method (linear-in-means)                       | HT    | Haj   | F     | L     | F- $\phi_0(G_1)$ | F- $\phi_0(G_2)$ | ND-F  | ND- $\phi_0(G_1)$ | ND- $G_1$ | ND-L  | ND- $\phi_0(G_2)$ | ND- $G_2$ |
|----------------------------------|-------|-------|-------|-------|------------------|------------------|-------|-------------------|-----------|-------|-------------------|-----------|
| **Empirical absolute bias**      | 0.001 | 0.001 | 0.000 | 0.000 | 0.001            | 0.001            | 0.000 | 0.001             | 0.895     | 0.000 | 0.002             | 0.998     |
| **Oracle SE**         | 0.087 | 0.086 | 0.068 | 0.068 | 0.055 | 0.054 | 0.068 | 0.055 | 0.217 | 0.068 | 0.056 | 0.676 |
| **Estimated SE**     | 0.086 | 0.085 | 0.067 | 0.067 | 0.054 | 0.054 | 0.067 | 0.053 | 0.053 | 0.067 | 0.054 | 0.053 |
| **Oracle coverage probability**  | 0.949 | 0.949 | 0.950 | 0.950 | 0.953            | 0.951            | 0.953 | 0.949             | 0.015     | 0.950 | 0.949             | 0.679     |
| **Empirical coverage probability**| 0.946 | 0.947 | 0.948 | 0.948 | 0.944            | 0.946            | 0.947 | 0.945             | 0.000     | 0.947 | 0.941             | 0.039     |









| Method (nonlinear contagion)     | HT    | Haj   | F     | L     | F- $\phi_0(G_1)$ | F- $\phi_0(G_2)$ | ND-F  | ND- $\phi_0(G_1)$ | ND- $G_1$ | ND-L  | ND- $\phi_0(G_2)$ | ND- $G_2$ |
|----------------------------------|-------|-------|-------|-------|------------------|------------------|-------|-------------------|-----------|-------|-------------------|-----------|
| **Empirical absolute bias**      | 0.000 | 0.000 | 0.000 | 0.000 | 0.000            | 0.000            | 0.000 | 0.000             | 0.176     | 0.000 | 0.000             | 0.044     |
| **Oracle SE**        | 0.052 | 0.020 | 0.018 | 0.018 | 0.017 | 0.017 | 0.018 | 0.017 | 0.073 | 0.018 | 0.017 | 0.212 |
| **Estimated SE**     | 0.052 | 0.021 | 0.020 | 0.019 | 0.018 | 0.018 | 0.020 | 0.018 | 0.018 | 0.019 | 0.018 | 0.017 |
| **Oracle coverage probability**  | 0.950 | 0.947 | 0.949 | 0.949 | 0.947            | 0.948            | 0.947 | 0.949             | 0.337     | 0.949 | 0.949             | 0.945     |
| **Empirical coverage probability**| 0.952 | 0.959 | 0.961 | 0.959 | 0.957            | 0.961            | 0.960 | 0.952             | 0.024     | 0.958 | 0.958             | 0.134     |

**Table 1**: Simulation results: network size $n=3000\left(b_n=3\right) . \tau=0.934$ for the linear-in-means setting and $\tau=0.193$ for the nonlinear contagion setting. We report the square of oracle SE and estimated SE in our tables for abbreviation.


## Real-world result:

| Method                          | HT    | Haj   | F     | L     | ND-F | ND- $\phi_0\left(G_1\right)$ | ND-L | ND- $\phi_0\left(G_2\right)$ |
|----------------------------------|-------|-------|-------|-------|------|----------------------------|-------|----------------------------|
| **Direct effect**                   | 0.0082 | 0.0146 | 0.0170 | 0.0168 | 0.0167 | 0.0161                     | 0.0164 | 0.0175                    |
| **Estimated SE**                 | 0.0272 | 0.0225 | 0.0218 | 0.0197 | 0.0211 | 0.0205                     | 0.0195 | 0.0193                    |
| **Spillover effect**                    | 0.0381 | 0.0611 | 0.0604 | 0.0581 | 0.0660 | 0.0686                     | 0.0592 | 0.0603                    |
| **Estimated SE**                 | 0.0447 | 0.0292 | 0.0270 | 0.0241 | 0.0258 | 0.0250                     | 0.0236 | 0.0233                    |

**Table 2**: Estimates $\hat{\tau}$ and Estimated SE of empirical experiments.












## Counter-example:
| Method                          | HT    | Haj   | F     | L     | ND-F | ND-L | ND- $\phi_0$ (G) |
|----------------------------------|-------|-------|-------|-------|------|------|-------------------|
| **Empirical absolute bias**      | 0.001 | 0.001 | 0.000 | 0.000 | 0.000 | 0.000 | 0.001             |
| **Oracle SE**                    | 0.037 | 0.037 | 0.041 | 0.050 | 0.037 | 0.043 | 0.038             |
| **Estimated SE**                 | 0.066 | 0.066 | 0.068 | 0.058 | 0.066 | 0.054 | 0.066             |
| **Oracle coverage probability**  | 0.950 | 0.950 | 0.950 | 0.950 | 0.950 | 0.950 | 0.950             |
| **Empirical coverage probability**| 0.999 | 0.999 | 0.999 | 0.978 | 0.999 | 0.986 | 0.999             |

**Table S.1**: The counter-example. $\tau=0.024$.






<!--

## Bandwidth selection and comparison:

We test the Bandwidth $=1,2,3,4$.

| Method (linear-in-means)         | HT            | Haj           | F             | L             | F- $\phi_0(G_1)$ | F- $\phi_0(G_2)$ | ND-F          | ND- $\phi_0(G_1)$ | ND- $G_1$     | ND-L          | ND- $\phi_0(G_2)$ | ND- $G_2$     |
|----------------------------------|----------------|----------------|----------------|----------------|------------------|------------------|----------------|-------------------|---------------|----------------|-------------------|---------------|
| **Empirical absolute bias**      | 0.935<br>0.934<br>0.934<br>0.935 | 0.935<br>0.933<br>0.933<br>0.934 | 0.935<br>0.933<br>0.933<br>0.934 | 0.935<br>0.933<br>0.933<br>0.934 | 0.935<br>0.934<br>0.933<br>0.933   | 0.935<br>0.935<br>0.934<br>0.933   | 0.935<br>0.933<br>0.933<br>0.934 | 0.935<br>0.934<br>0.933<br>0.933     | 0.0463<br>0.044<br>0.037<br>0.039 | 0.935<br>0.933<br>0.933<br>0.934 | 0.935<br>0.935<br>0.934<br>0.934     | -0.072<br>-0.055<br>-0.067<br>-0.070 |
| **Oracle SE**                    | 0.085<br>0.085<br>0.087<br>0.085 | 0.084 <br>0.084<br>0.086<br>0.084 | 0.067<br>0.066<br>0.068<br>0.067 | 0.067<br>0.066<br>0.068<br>0.067 | 0.055<br>0.055<br>0.055<br>0.054   | 0.054<br>0.054<br>0.055<br>0.053   | 0.067<br>0.066<br>0.068<br>0.067 | 0.055<br>0.055<br>0.056<br>0.054     | 0.021<br>0.216<br>0.217<br>0.222 | 0.067<br>0.066<br>0.068<br>0.067 | 0.055<br>0.055<br>0.055<br>0.054     | 13.46<br>0.666<br>0.674<br>0.680 |
| **Estimated SE**                 | 0.086<br>0.087<br>0.087<br>0.087 | 0.084<br>0.085<br>0.085<br>0.085 | 0.067<br>0.067<br>0.067<br>0.067| 0.067<br>0.067<br>0.067<br>0.067 | 0.054<br>0.054<br>0.054<br>0.054   | 0.054<br>0.054<br>0.054<br>0.054   | 0.067<br>0.067<br>0.067<br>0.067 | 0.054<br>0.054<br>0.054<br>0.054     | 0.053<br>0.053<br>0.053<br>0.053 | 0.067<br>0.067<br>0.067<br>0.067 | 0.054<br>0.053<br>0.053<br>0.053     | 0.053<br>0.053<br>0.053<br>0.053 |
| **Oracle coverage probability**  | 0.951<br>0.951<br>0.950<br>0.950 | 0.951<br>0.951<br>0.951<br>0.952 | 0.951<br>0.950<br>0.949<br>0.950 | 0.950<br>0.950<br>0.949<br>0.951 | 0.951<br>0.950<br>0.949<br>0.951   | 0.950<br>0.951<br>0.948<br>0.951   | 0.950<br>0.952<br>0.949<br>0.951 | 0.951<br>0.951<br>0.950     | 0.012<br>0.016<br>0.015<br>0.019 | 0.950<br>0.951<br>0.948<br>0.951 | 0.950<br>0.951<br>0.951<br>0.950     | 0.991<br>0.681<br>0.680<br>0.681 |
| **Empirical coverage probability**| 0.953<br>0.952<br>0.949<br>0.954 | 0.953<br>0.951<br>0.948<br>0.954 | 0.951<br>0.953<br>0.947<br>0.952 | 0.950<br>0.953<br>0.947<br>0.952 | 0.946<br>0.945<br>0.941<br>0.948   | 0.948<br>0.949<br>0.944<br>0.951   | 0.950<br>0.953<br>0.946<br>0.951 | 0.945<br>0.944<br>0.940<br>0.949     | 0.000<br>0.000<br>0.000<br>0.000 | 0.950<br>0.953<br>0.946<br>0.950 | 0.944<br>0.942<br>0.939<br>0.947     | NA<br>0.041<br>0.039<br>0.042 |




-->














#### Additional visualisation:



| Network | Covariate X | Noise ε |
|-------------|---------|-----------|
| ![](figures/output(1).png) | ![](figures/output(2).png) |![](figures/output(3).png) |


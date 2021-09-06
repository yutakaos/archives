# siarbeta : Beta-dependent Stable Isotope Mixing Model

A tutorial for siarbeta package.

## Installation
```r
library(devtools)
devtools::install_github("yutakaos/archives/simm/siarbeta")
```

## Tutorial

### Load library and set data
```r
# Load libarary
library(siarbeta)

# Set data
library(siar)
data("geese1demo", "sourcesdemo", "correctionsdemo", "concdepdemo")
mixture <- geese1demo
sources <- sourcesdemo[,-1]
correct <- correctionsdemo[,-1]
concdep <- concdepdemo[,-1]
source_names <- as.character(sourcesdemo[,1])
rm(geese1demo, sourcesdemo, correctionsdemo, concdepdemo)
```

### Run MCMC
```r
# Ordinary SIAR
out_L = siarbeta(
    mixture, sources, correct, concdep, alpha = 1, beta = 1,
    error = "parnell", source_names = source_names,
    chains = 3, iters = 5000, burns = 1000, thins = 4)

# Beta-dependent SIAR
out_H = siarbeta(
    mixture, sources, correct, concdep, alpha = 1, beta = 1000,
    error = "parnell", source_names = source_names,
    chains = 3, iters = 5000, burns = 1000, thins = 4)

# Summary statistics
summary(out_L)
#                      Mean         SD          2.5%          50%       97.5%     rhat   mean_r2
# Zostera        0.59970605 0.10624537   0.399537936   0.59719530   0.8102002 1.000347 0.1858890
# Grass          0.07019535 0.02377007   0.024779325   0.06929327   0.1166252 1.001158 0.1215366
# U.lactuca      0.13130782 0.09348966   0.005846707   0.11301847   0.3465096 1.001306 0.1330593
# Enteromorpha   0.19879078 0.13265536   0.008563942   0.18469323   0.4851537 1.000873 0.3301500
# SD1            0.64726213 0.39427916   0.104439547   0.57393721   1.6030457 1.000000        NA
# SD2            1.02409361 0.54679174   0.218837011   0.94059419   2.2720622 1.001462        NA
# lp           -26.10691496 1.66110829 -30.268261223 -25.78067618 -23.9629372 1.000182        NA

summary(out_H)
#                      Mean          SD          2.5%          50%        97.5%     rhat   mean_r2
# Zostera        0.57560782 0.008385274   0.559376168   0.57550796   0.59261019 1.009336 0.6277728
# Grass          0.05774041 0.002534191   0.053211424   0.05767994   0.06306487 1.004999 0.7392251
# U.lactuca      0.02632429 0.015375535   0.001734123   0.02535515   0.05942069 1.005888 0.7282766
# Enteromorpha   0.34032748 0.024441781   0.288189375   0.34121254   0.38213397 1.007661 0.8440318
# SD1            0.01791682 0.009200850   0.003226429   0.01684612   0.03767967 1.000634        NA
# SD2            0.68406997 0.017655317   0.647955216   0.68419640   0.71707637 1.000356        NA
# lp           -23.50743517 0.008424708 -23.525708883 -23.50549092 -23.50184471 1.037301        NA

# Correlation plot
plot_corr(out_L[,1:4])
```
<figure>
<img src="tools/figures/Fig_4.png" width="70%">
<figcaption><i>Figure 1 | Correlation plot.</i></figcaption>
</figure>

### Diagnose underdetermination
```r
# Diagnose underdetermination
evaluate_ump(out_H)
#                      lower        upper       delta
# Zostera       5.585727e-01   0.59172433 0.033151640
# Grass         5.305615e-02   0.06271182 0.009655665
# U.lactuca     1.674532e-05   0.05427174 0.054254994
# Enteromorpha  2.909324e-01   0.38310595 0.092173571
# SD1           1.929180e-03   0.03496760 0.033038423
# SD2           6.474285e-01   0.71615196 0.068723415
# lp           -2.351821e+01 -23.50117919          NA

# Plot mixing isotope space
mixingspace(
    mixture, sources, correct, axis = 2:1,
    source_names  = source_names,
    element_names = c("d13C", "d15N") )

# Posterior density plot
plot_post(out_L, out_H, type = "source")
```
<figure>
<img src="tools/figures/Fig_3.png" width="70%">
<figcaption><i>Figure 2 | Isotopic mixing space and posterior density plot.</i></figcaption>
</figure>


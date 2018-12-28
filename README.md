Selecting representative subsets of vignettes for investigating multiple facets of moral judgement
==================================================================================================
This repository contains the Matlab code used by Kruepke et al. (2018) to generate subsets of 39 moral vignettes from the set of 312 moral vignettes created by Knutson et al. (2010). Data from Knutson et al. (2010) gives the 312 vignettes with ratings on 13 features (e.g., emotional intensity, emotional aversion, harm, self-benefit, etc.) by 30 subsets. The examples below use the Matlab code in this repository to produce surveys with 39 vignettes by selecting one vignette that is labeled "low", "moderate/neutral", and "high" for each of the 13 features. 


Example 1
---------
In this example, vignettes are labeled "low", "moderate/neurtral", and "high" for each of the 13 features based on whether the mean rating across subjects for a given feature falls within the intervals [1, 3], (3,4), and [4, 7], respectively.
```
seed = 12345;

nlow = 1; nmid = 1; nhigh = 1;
breakdown = [nlow, nmid, nhigh];

a = 1.0; b = 3.0; d = 4.0; e = 7.0;
limits = [a, b, d, e];

select_vignette_subset('data/SuppData.csv', 'subQ_limits', limits, breakdown, 1, 'heuristic', seed);
```


Example 2
---------
In this example, vignettes are labeled "low", "moderate/neurtral", and "high" for each of the 13 features based on whether the mean rating across subjects for a given feature falls within the respective intervals: [min, u-s], (u-s, u+s), [u+s, max], where min is the minimum rating for a given feature across all 312 vignettes, u is the mean, s is the standard deviation, and max is the maximum.
```
seed = 12345;

nlow = 1; nmid = 1; nhigh = 1;
breakdown = [nlow, nmid, nhigh];

select_vignette_subset('knuston_2010_data.csv', 'subQ', [], breakdown, 3, 'heuristic', seed);
select_vignette_subset('knuston_2010_data.csv', 'subQ', [], breakdown, 3, 'random', seed);
```


Citations
---------
Knutson et al. (2010). Behavioral norms for condensed moral vignettes. *Social Cognitive and Affective Neuroscience* 5(4): 378-384. [doi:10.1093/scan/nsq005](https://doi.org/10.1093/scan/nsq005)

Kruepke, M. D., Molloy, E. K., Bresin, K. W., Barbey, A. K., and Verona, E. (2017). A Brief Assessment Tool for Investigating Facets of Moral Judgment from Realistic Vignettes. *Behavior Research Methods*. 50(3):922-936. [doi:10.3758/s13428-017-0917-3](https://doi.org/10.3758/s13428-017-0917-3)

Molloy, E. K. and Kruepke, M. D., (2017). Selecting representative subsets of vignettes for investigating multiple facets of moral judgement: Documentation and MATLAB Code. GitHub repository, [https://github.com/ekmolloy/vignette-selection](https://github.com/ekmolloy/select-vignette-subsets).
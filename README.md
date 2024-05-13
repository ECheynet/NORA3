# Gridded NORA3 data: automated and remote data extraction 

[![View Gridded NORA3 data: automated and remote data extraction on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/93685-gridded-nora3-data-automated-and-remote-data-extraction)

[![Donation](https://camo.githubusercontent.com/a37ab2f2f19af23730565736fb8621eea275aad02f649c8f96959f78388edf45/68747470733a2f2f77617265686f7573652d63616d6f2e636d68312e707366686f737465642e6f72672f316339333962613132323739393662383762623033636630323963313438323165616239616439312f3638373437343730373333613266326636393664363732653733363836393635366336343733326536393666326636323631363436373635326634343666366536313734363532643432373537393235333233303664363532353332333036313235333233303633366636363636363536353264373936353663366336663737363737323635363536653265373337363637)](https://www.buymeacoffee.com/echeynet)

## Summary
The NORA3 atmosphere hindcast data are remotely accessed using the OPeNDAP framework and the THREDDS Data Server of the Norwegian Meteorological Institute [1]. The data can be extracted for multiple latitudes and longitudes and stored in a gridded format. Atmospheric data are retrieved for seven different heights, from 10 m to 750 m above the surface. The mean wind speed profile can be interpolated using a non-linear scheme based on the Deaves and Harris model [2,3]. Above the ocean, sea roughness is modeled using the Charnock relation. A summary of the NORA3 data and their potential is available in [4]. This code was used to extract the data used in [5].

The `getNORA3_subset` function specifically accesses subsets of the NORA3 data (monthly storage), increasing data collection by a factor of x10 or more. However, these subsets do not contain profiles of temperature.

## Content

The present submission contains:
  - `getNORA3.m`: A function that imports the NORA3 hindcast from [1].
  - `getNORA3_subset.m`: A function accessing subsets of NORA3 for enhanced data collection. Note: This subset does not include temperature profiles.
  - `getz0_charnock.m`: A function estimating sea surface roughness using Charnock's equation and the logarithmic profile of the mean wind speed.
  - `interpU.m`: A function interpolating the mean wind speed using The Deaves and Harris model [2,3].
  - `world.mat`: A data file containing the coastline of countries (used for visualization purposes only).
  - `Documentation.mlx`: A Matlab LiveScript illustrating how these functions can be used.
  - `Documentation2.mlx`: A second Matlab LiveScript illustrating how these functions can be used.
  - `Documentation3-subset.mlx`: A third Matlab LiveScript illustrating getNORA3_subset works. Note that the surface temperature seems a little different here.

Any comments, questions, or suggestions are welcome.

## References

[1] https://thredds.met.no  

[2] Harris, R. I., & Deaves, D. M. (1980, November). The structure of strong winds, wind engineering in the eighties. In Proc. CIRIA Conf.

[3] ESDU. (1985). ESDU 85020-Characteristics of atmospheric turbulence near the ground. Part II: single point data for strong winds (neutral atmosphere).

[4] Solbrekke, I. M., Sorteberg, A., & Haakenstad, H. (2021). Norwegian hindcast archive (NORA3)–A validation of offshore wind resources in the North Sea and Norwegian Sea. Wind Energy Science Discussions, 1-31.

[5] Cheynet, E., Solbrekke, I. M., Diezel, J. M., & Reuder, J. A one-year comparison of new wind atlases over the North Sea. Journal of Physics: Conference Series 2362 (1), 012009


## Example

<img src="illustration.jpg" alt="Mean wind speed at 10 m above the surface in Northern Europe" width="700"/>

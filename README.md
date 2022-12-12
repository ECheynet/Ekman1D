# Ekman1D
Variable eddy viscosity Ekman layer in the ABL (1D)

[![View Variable eddy viscosity Ekman layer in the ABL (1D) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/74316-variable-eddy-viscosity-ekman-layer-in-the-abl-1d)
[![DOI](https://zenodo.org/badge/264255959.svg)](https://zenodo.org/badge/latestdoi/264255959)


[![Donation](https://camo.githubusercontent.com/a37ab2f2f19af23730565736fb8621eea275aad02f649c8f96959f78388edf45/68747470733a2f2f77617265686f7573652d63616d6f2e636d68312e707366686f737465642e6f72672f316339333962613132323739393662383762623033636630323963313438323165616239616439312f3638373437343730373333613266326636393664363732653733363836393635366336343733326536393666326636323631363436373635326634343666366536313734363532643432373537393235333233303664363532353332333036313235333233303633366636363636363536353264373936353663366336663737363737323635363536653265373337363637)](https://www.buymeacoffee.com/echeynet)


## Summary

Matlab implementation of the solution to the Ekman equations in the atmospheric boundary layer. The flow is assumed horizontal and homogeneous. however, a height-dependant eddy viscosity can be modelled. The solutions are provided in one-dimension.

## Content

The submission includes

  - The function EkmanAnalytic that provides analytics solution of Ekman's equations for a constant eddy viscosity in the atmospheric boundary layer.
  - The function solveEkman that numerically solves Ekman's equations with an explicit finite difference scheme and allows the use of height-dependant eddy viscosity. The numerical implementation is partly inspired by [1].
  - An example file Example0.mlx and reproduces some of the figures displayed in ref [2]

Any question, suggestion or comment is welcomed.

## References:

[1] Berger, B. W., & Grisogono, B. (1998). The baroclinic, variable eddy viscosity Ekman layer. Boundary-layer meteorology, 87(3), 363-380.

[2] https://www.whoi.edu/science/PO/people/jprice/website/education_scripts.html

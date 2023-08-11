# Ekman1D
Variable eddy viscosity Ekman layer in the ABL (1D)

[![View Variable eddy viscosity Ekman layer in the ABL (1D) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/74316-variable-eddy-viscosity-ekman-layer-in-the-abl-1d)
[![DOI](https://zenodo.org/badge/264255959.svg)](https://zenodo.org/badge/latestdoi/264255959)

<a href="https://www.buymeacoffee.com/echeynet" target="_blank"><img src="https://www.buymeacoffee.com/assets/img/custom_images/orange_img.png" alt="Buy Me A Coffee" style="height: 25px !important;width: 120px !important;box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;-webkit-box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;" ></a>


## Summary

Matlab implementation of the solution of the Ekman equations in the atmospheric boundary layer. The flow is assumed to be horizontal and homogeneous, but a height dependent turbulent viscosity can be modelled. The solutions are provided in one dimension.

## Contents

The submission contains

  - The function EkmanAnalytic, which provides an analytical solution to Ekman's equations for a constant eddy viscosity in the atmospheric boundary layer.
  - The solveEkman function, which numerically solves Ekman's equations using an explicit finite difference scheme and allows the use of height-dependent eddy viscosity. The numerical implementation is partly inspired by [1].
  - The function solveEkman_bcp4v.m, which uses the Matlab function bcp4v to solve Ekman's equations.
  - The function scm_bcp4v.m, which is used for more advanced single column models.
  - An example file Example0.mlx which reproduces some of the figures shown in ref [2].
  - An example file documentation_part2.mlx.mlx which uses the scm_bcp4v function.

Any questions, suggestions or comments are welcome.

## References:

[1] Berger, B. W., & Grisogono, B. (1998). The baroclinic, variable eddy viscosity Ekman layer. Boundary-layer meteorology, 87(3), 363-380.

[2] https://www.whoi.edu/science/PO/people/jprice/website/education_scripts.html

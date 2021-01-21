# MiMA [![DOI](https://zenodo.org/badge/36012278.svg)](https://zenodo.org/badge/latestdoi/36012278)
Model of an idealized Moist Atmosphere

MiMA is an intermediate-complexity General Circulation Model with interactive water vapor and full radiation. It is publicly available but users are asked to cite the appropriate references in any publication resulting from the use of MiMA. Please refer to the below sections for more information.

* [Getting started](GettingStarted.md): Downloading, compiling, test run
* [Parameter settings](Parameters.md): Default values and settings
* [Version history](Versions.md): History of main additions and changes
* [References](#references): Required and relevant references
* [License](#license): License information

See the 30 second trailer on [YouTube](https://www.youtube.com/watch?v=8UfaFnGtCrk "Model of an idealized Moist Atmosphere (MiMA)"):

[![MiMA thumbnail](https://img.youtube.com/vi/8UfaFnGtCrk/0.jpg)](https://www.youtube.com/watch?v=8UfaFnGtCrk "Model of an idealized Moist Atmosphere (MiMA)")

## News
With time, I hope to post usage, papers, etc. which relate to MiMA here. If you're a MiMA user, please let me know of any news at <coding@martinjucker.com>.

* February 2020: New release (v1.1) adds restart from arbitrary initial conditions (as used in [Yamada and Pauluis (2017)](https://doi.org/10.1175/JAS-D-16-0329.1) below) and possibility to use the Navy high-resolution land-sea mask.
* January 2020: Two new papers published using MiMA for research on climate dynamics. These papers have further developed MiMA and their additions will be published in the upcoming v2.0: 
  * [Garfinkel, C.I., I. White, E.P. Gerber, M. Jucker, and M. Erez (2020): _The building blocks of Northern Hemisphere wintertime stationary waves_, Journal of Climate, doi:10.1175/JCLI-D-19-0181.1](https://doi.org/10.1175/JCLI-D-19-0181.1) which, among other important findings, shows that MiMA is just as good as any CMIP5 model to represent Northern Hemisphere stationary waves.
  * [White, I.P., C. Garfinkel, E. Gerber, M. Jucker, P. Hitchcock, and J. Rao (2020): _The generic nature of the tropospheric response to sudden stratospheric warmings_, Journal of Climate, doi: 10.1175/JCLI-D-16-0703.1](https://doi.org/10.1175/JCLI-D-16-0703.1) which examines the evolution of naturally occurring versus synthetically forced Sudden Stratospheric Warmings.
* August 2019: [Iceberg](https://research-iceberg.github.io) article about MiMA's various choices for mixed layer depth and albedo distribution: [_The surface of an aquaplanet GCM_](https://research-iceberg.github.io/papers/M_Jucker_201907/). PDF: <a href="https://doi.org/10.5281/zenodo.3358284" target="_blank">doi.org/10.5281/zenodo.3358284</a>
* March 2018: [Isca framework](https://empslocal.ex.ac.uk/people/staff/gv219/ISCA/index.html) (which contains MiMA) reference paper published in Geoscientific Model Development: [Vallis _et al._, 2018: *Isca, v1.0: a framework for the global modelling of the atmospheres of Earth and other planets at varying levels of complexity*, Geosci. Model Dev., doi:10.5194/gmd-11-843-2018](https://doi.org/10.5194/gmd-11-843-2018)
* February 2018: Talk and poster about MiMA at the [2nd Pan-GASS meeting](http://singh.sci.monash.edu/Pan-GASS/index.shtml). The poster can be viewed [here](https://doi.org/10.6084/m9.figshare.6118934.v1).
* August 2017: Bugfix patch (v1.0.1) addresses incoming SW radation bugs.
* June 2017: Ray Yamada and Olivier Pauluis use MiMA for baroclinic lifecycle studies: [Yamada, R., and O. Pauluis, 2017: Wave-mean-flow interactions in moist baroclinic lifecycles. J. Atmo. Sci., 74, 2143-2162, doi:10.1175/JAS-D-16-0329.1](https://doi.org/10.1175/JAS-D-16-0329.1)
* June 2017: MiMA reference paper published in Journal of Climate:
[M Jucker and EP Gerber, 2017: *Untangling the annual cycle of the tropical tropopause layer with an idealized moist model*, J Clim, doi:10.1175/JCLI-D-17-0127.1](http://dx.doi.org/10.1175/JCLI-D-17-0127.1).



## References

MiMA
* [Jucker and Gerber, J Clim (2017)](http://dx.doi.org/10.1175/JCLI-D-17-0127.1)
* [Garfinkel et al., J Clim (2020)](http://journals.ametsoc.org/doi/10.1175/JCLI-D-19-0181.1)
* code:
  * latest version: [![DOI](https://zenodo.org/badge/36012278.svg)](https://zenodo.org/badge/latestdoi/36012278)
  * v1.1: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3637607.svg)](https://doi.org/10.5281/zenodo.3637607)
  * v1.0: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.321708.svg)](https://doi.org/10.5281/zenodo.321708)

Gray radiation model
* [Frierson, Held, Zurita-Gotor, JAS (2006)](http://journals.ametsoc.org/doi/abs/10.1175/JAS3753.1)
* [Frierson, JAS (2007)](http://journals.ametsoc.org/doi/abs/10.1175/JAS3935.1)
* [Frierson, Held, Zurita-Gotor, JAS (2007)](http://journals.ametsoc.org/doi/abs/10.1175/JAS3913.1)

RRTM
* [Mlawer et al., JGR (1997)](http://doi.wiley.com/10.1029/97JD00237)
* [Iacono et al., JGR (2000)](http://doi.wiley.com/10.1029/2000JD900091)
* [Iacono et al., JGR (2008)](http://onlinelibrary.wiley.com/doi/10.1029/2008JD009944/abstract)
* [Clough et al., JQSRT (2005)](http://www.sciencedirect.com/science/article/pii/S0022407304002158)


## License

MiMA is distributed under a GNU GPLv3 license. That means you have permission to use, modify, and distribute the code, even for commercial use. However, you must make your code publicly available under the same license. See LICENSE.txt for more details.

AM2 is distributed under a GNU GPLv2 license. That means you have permission to use, modify, and distribute the code, even for commercial use. However, you must make your code publicly available under the same license.

RRTM/RRTMG: Copyright © 2002-2010, Atmospheric and Environmental Research, Inc. (AER, Inc.). This software may be used, copied, or redistributed as long as it is not sold and this copyright notice is reproduced on each copy made. This model is provided as is without any express or implied warranties.

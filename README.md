## DRIFTER SIMULATION TOOLBOX FOR MATLAB ##

v9.0

[![DOI](https://zenodo.org/badge/198690464.svg)](https://zenodo.org/badge/latestdoi/198690464)

Tools for simulating drifters given gridded vector data, such as from 
oceanographic HF radar data. Use of this code has resulted in at least two papers, which are cited below.
Improvements developed for Matson et al. (2019), which includes error estimation, are not yet included here. 
Otherwise the toolbox includes:   

- Use of 4th order Runge-Kutta for accurate integration of discrete velocities.
- Optional settings to simulate subgrid turbulence.
- Dependencies in the /private folder. Also see the note about the use of code folding and subfunctions.

HOW TO USE IT
- download it and add the unzipped directory to your MATLAB path
- start with the forthcomming example ...


AUTHORS  
Brian Emery  
Krisada Letcherononguong  
Kirk Ireson  
Chris Gotschalk  
Carter Ohlmann  
Libe Washburn  


ACKNOWLEDGMENT

The development of this software was made possible in part through support from the National Science
Foundation under Grant OCE-0352187 and the Minerals Management Service, U.S. Department of Interior
under MMS Agreement 1435-01-00-CA-31063 to 18212. The views and conclusions contained herein are
those of the authors and should not be interpreted as necessarily representing the official 
policies, either expressed or implied, of the U.S. government.


REFERENCES

Ohlmann, J. C., J. H. LaCasce, L. Washburn, A. J. Mariano, and B. Emery, 2012,
Relative dispersion observations and trajectory modeling in the 
Santa Barbara Channel, J. Geophys. Res., 117, C05040, doi:10.1029/2011JC007810.

Matson, Paul G., L. Washburn, E. A. Fields, C. Gotschalk, T. M. Ladd, D. A. 
Siegel, Z. S. Welch, and M. D. Iglesias‐Rodriguez. "Formation, development, 
and propagation of a rare coastal coccolithophore bloom." Journal of 
Geophysical Research: Oceans 124, no. 5 (2019): 3298-3316.


NOTES

- .m files use functions as blocks of code, so code folding (for example using the 
  keyboard shortcut 'cmd =' ) makes it easy to move among these.
- Data structures contain variables of similar origin following the HFRProgs
  convention, rows = locations, cols = time
  
HOW TO CITE
  
```
@misc{Emery2018code,
  author       = {Emery, Brian and Letcherononguong, Krisada and Ireson, Kirk and Gotschalk, Chris and Ohlmann, Carter and Washburn, Libe},
  title        = {{Drifter Simulation Toolbox for MATLAB}, Software Release Version 9.0, https://doi.org/10.5281/zenodo.3350830},
  version      = {9.0},
  month        = Jul,
  year         = 2019,
  doi          = {10.5281/zenodo.3350830},
  url          = {https://doi.org/10.5281/zenodo.3350830}
}
```




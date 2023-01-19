
[![Build Status](https://dev.azure.com/mdolab/Private/_apis/build/status/mdolab.geograd?repoName=mdolab%2Fgeograd&branchName=master)](https://dev.azure.com/mdolab/Private/_build/latest?definitionId=21&repoName=mdolab%2Fgeograd&branchName=master)

# GeoGrad - A spatial integration constraint for MDO 
GeoGrad is a general geometric constraint for spatial integration in gradient-based optimization.
It can be used for aerodynamic shape optimization within the [MACH framework](https://github.com/mdolab/MACH-Aero) through [pyGeo](https://github.com/mdolab/pygeo).

## Documentation
See the details in the `report` directory for information on the parallelism implemented in GeoGrad. 

## Citation
Please cite GeoGrad by referencing [this paper](http://www.umich.edu/~mdolaboratory/pdf/Brelje2020a.pdf).

B. J. Brelje, J. Anibal, A. Yildirim, C. A. Mader, and J. R. R. A. Martins. AIAA Journal, 58(6):2571–2580, 2020.
```
@article{Brelje2020a,
  author = {Brelje, Benjamin J. and Anibal, Joshua and Yildirim, Anil and Mader, Charles A. and Martins, Joaquim R. R. A.},
  doi = {10.2514/1.J058366},
  journal = {AIAA Journal},
  keywords = {ank, ccavd},
  month = jun,
  number = {6},
  pages = {2571--2580},
  title = {Flexible Formulation of Spatial Integration Constraints in Aerodynamic Shape Optimization},
  volume = {58},
  year = {2020}
}
```

## License
GeoGrad is licensed under the Apache License, Version 2.0 (the "License"). See `LICENSE` for the full license.


### Latest
- Vectorized AABB tests to exclude swathes of computation
- Tapenade for the main routines
- Dynamic load balancing
- Need to invent an optimization case where the load balancing actually changes

# AtmosFOAM-multiFluid

Tools for modelling the atmosphere and atmospheric convection with multi-fluids. To be used with [AtmosFOAM-tools](https://github.com/AtmosFOAM/AtmosFOAM-tools/) and [AtmosFOAM](https://github.com/AtmosFOAM/AtmosFOAM/)


* Install [OpenFOAM dev](https://github.com/OpenFOAM/OpenFOAM-dev).
* Ensure [AtmosFOAM-tools](https://github.com/AtmosFOAM/AtmosFOAM-tools/) is installed
* Ensure [AtmosFOAM](https://github.com/AtmosFOAM/AtmosFOAM/) is installed
* Compile all AtmosFOAM-multiFluid applications and libraries using `./Allwmake`
* Export environment variables  in your `~/.bashrc` file:

       export ATMOSFOAM_TOOLS_SRC=/path/to/AtmosFOAM-tools/src
       export GMTU=/path/to/AtmosFOAM-tools/gmtUser
       export ATMOSFOAM_SRC=/path/to/AtmosFOAM/src
       export ATMOSFOAM_MULTIF_SRC=/path/to/AtmosFOAM-multiFluid/src


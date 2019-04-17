# rheoConvective
# FiniteVolume  infinite-Pr  Rayleigh-Bénard Solver
---------------------

OpenFOAM code which integrates the Rayleigh-Bénard convection equations in the Boussinesq approximation.
Infinite Pr number and Pressure&Temperature viscosity dependence is introduced. 

---------------------

Main equations 

div(U) = 0              ---- (PISO algorithm)

ddt(U) + div( U x U ) = -grad(p) + div(stress) - beta(T-TRef) g

ddt(T) + div(U T) = div(diffusion)


where beta is the thermal expansion.    


Dynamic (newtonian) viscosity 

		eta = 1/2A * exp((E+pV)/RT)


The test case is a 3D version of the 2D one described 
in Fowler et al. 2015; https://doi.org/10.1063/1.4923061


How to use convectiveFoam?

- clone the source files in your preferred directory:
    git clone https://github.com/user/rheoConvectiveFoam.git

- Be sure the environment variables of OpenFoam 5.x have been
  loaded in your terminal window

- enter in the rheoConvectiveFoam folder and install it running 
  the script:
  
    ./Allwmake

- try it with the Case, in your tutorial folder, by
  using the following commands:
  
    cd $WM_PROJECT_USER_DIR/tutorials/rheoConvectiveFoam/Case
    
    ./Allrun

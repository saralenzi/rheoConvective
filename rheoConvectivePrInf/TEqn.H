{
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
     ==
        fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();
    

    /////Add here------------------------------------
    //nu = nu0*exp((E+p*V)/(Rgas*T));
    //nu = nu0;
    ////---------------------------------------------
                            

    fvOptions.correct(T);
}

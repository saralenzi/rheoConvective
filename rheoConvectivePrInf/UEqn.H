//update viscosity 

//nu = nu0*exp((E+rho0*( mag(g)*z0 + (g & mesh.C()) )*Vol)/(Rgas*(T)));
nu = nu0*exp((E+rho0*( mag(g)*z0 - (g & mesh.C()) )*Vol)/(Rgas*(T)));
Info << "Max and min eta: " << max(nu).value() << "\t" << min(nu).value() << endl;


// Solve the Momentum equation
MRF.correctBoundaryVelocity(U);

dimensionedScalar Dt("Dt", dimTime, runTime.deltaTValue());
tmp<fvVectorMatrix> tUEqn
(
//    fvm::ddt(U) + fvm::div(phi, U)
//  + MRF.DDt(U)+ 
//  //turbulence->divDevReff(U)
//    turbulence->nu()*fvm::laplacian(U)

//    (fvm::Sp(scalar(1)/Dt, U) - U/Dt)
  - fvm::laplacian(nu, U)
  - fvc::div(nu*dev(fvc::grad(U)().T()))
 ==
    fvOptions(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    //for(int ii=0; ii<10; ii++)
    {
        solve
        (
            UEqn
         ==
           -fvc::grad(p)
          + fvc::reconstruct
            (
                -fvc::interpolate(beta*(T-TRef))*(g & mesh.Sf())
            )
        );
    }

    fvOptions.correct(U);
}

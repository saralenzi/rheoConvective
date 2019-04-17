#include "readGravitationalAcceleration.H"

Info<< "Reading thermophysical properties\n" << endl;

IOdictionary transportProperties
(
     IOobject
         (
	  "transportProperties",
	   runTime.constant(),
	   mesh,
	   IOobject::MUST_READ_IF_MODIFIED,
	   IOobject::NO_WRITE
          )
);



dimensionedScalar Rgas
(
      "Rgas",
       transportProperties.lookup("Rgas")
);

dimensionedScalar nu0
(
       "nu0",
        transportProperties.lookup("nu0")
);

dimensionedScalar E
(
       "E",
        transportProperties.lookup("E")
);


dimensionedScalar Vol
(
       "Vol",
        transportProperties.lookup("Vol")
);

dimensionedScalar rho0
(
       "rho0",
        transportProperties.lookup("rho0")
);

dimensionedScalar z0
(
       "z0",
        transportProperties.lookup("z0")
);




Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);




Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


    Info<< "E=" << E<<endl;
    Info<< "Rgas=" << Rgas<<endl;
    Info<< "nu0=" << nu0<<endl;
    Info<< "T=" << T<<endl;


Info<< "Reading field nu\n" << endl;
volScalarField nu
(
    IOobject
    (
       "nu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    nu0*exp((E+rho0*(mag(g)*z0 + (g & mesh.C()))*Vol)/(Rgas*T))
);
nu.write();


#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

// Thermal expansion coefficient [1/K]
dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    laminarTransport
);

// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);

// Laminar Prandtl number
dimensionedScalar Pr("Pr", dimless, laminarTransport);

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, laminarTransport);

Info<< "Creating turbulence model\n" << endl;
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// kinematic turbulent thermal thermal conductivity m2/s
Info<< "Reading field alphat\n" << endl;
volScalarField alphat
(
    IOobject
    (
        "alphat",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);




#include "createMRF.H"





Info<< "Reading transportProperties\n" << endl;

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

dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties
);

dimensionedScalar mu
(
    "mu",
    dimensionSet(1, 1, -2, 0, 0, -2, 0),
    transportProperties
);

dimensionedScalar sigma
(
    "sigma",
    dimensionSet(-1, -3, 3, 0, 0, 2, 0),
    transportProperties
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
/*
Info<< "Deffining Si\n" << endl;
volScalarField Si
(
    IOobject
    (
        "Si", //dimensionSet(0, 0, 0, 0, 0, 0, 0),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
); */

Info<< "Reading field rho\n" << endl; //I added the temporal dependence for rho
volScalarField rho
(
     IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh, //thermo.rho()
    dimensionedScalar("rho",dimensionSet(1,-3,0,0,0,0,0),10)
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

#include "createPhi.H"

Info<< "Reading field pB\n" << endl;
volScalarField pB
(
    IOobject
    (
        "pB",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field B\n" << endl;
volVectorField B
(
    IOobject
    (
        "B",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


#include "createPhiB.H"

dimensionedScalar DB = 1.0/(mu*sigma);
DB.name() = "DB";

dimensionedScalar DBU = 1.0/(2.0*mu); //I removed rho from divisor
DBU.name() = "DBU";


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, piso.dict(), pRefCell, pRefValue);

mesh.setFluxRequired(p.name());
mesh.setFluxRequired(pB.name());

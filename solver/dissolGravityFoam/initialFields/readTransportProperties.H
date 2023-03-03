Info<< "==Reading  transport Properties==\n" << endl;

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

dimensionedScalar Df
(
    "Df",
    dimViscosity,
    transportProperties.lookup("Df")
);
Info<< "effective molecular diffusivity: "<<Df<<endl;

dimensionedScalar rho0
(
    "rho0",
    dimDensity,//
    transportProperties.lookup("rho0")
);

dimensionedScalar beta0
(
    "beta0",
	dimensionSet(-1, 3, 0, 0, 0, 0, 0),
    //dimDensity/(dimDensity*dimDensity),//dimless/dimDensity,
    transportProperties.lookup("beta0")
);

dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties.lookup("nu")
);

Info<< "==Complete reading transport Properties== "<<nl<<endl;
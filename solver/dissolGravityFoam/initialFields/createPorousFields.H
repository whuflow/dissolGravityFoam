
Info<< "==Reading porous medium transportProperties==\n" << endl;

IOdictionary porousTransportProperties
(
    IOobject
    (
        "porousTransportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar K0
{
    "K0",
    dimensionSet(0, 2, 0, 0, 0, 0, 0),
    porousTransportProperties.lookup("K0")
};
Info<< "permeability property K0: "<<K0<<endl;

dimensionedScalar rK0
{
    "rK0",
    dimensionSet(0, -2, 0, 0, 0, 0, 0),
    (1.0/K0).value()
};

Info<< "Reading field porosity" << endl;
volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

scalar minEps=1e-4;
volScalarField tEps(max(eps,minEps));
volScalarField rEps
{
    "rEps",
    1.0/tEps
};

surfaceScalarField rEpsf
(
    "repsf",
    fvc::interpolate(1/tEps)
);

Info<< "creating solid fraction" << endl;

volScalarField epsSolid
(
    "epsSolid",
    1-eps
);

Info<< "create permeability field" << endl;


volScalarField rK
(
    "rK",
    rK0*(1.0-eps)*(1.0-eps)/(tEps*tEps*tEps)
);
// rK.write();

surfaceScalarField nEpsv
(
   IOobject
   (
      "nEpsf",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
   ),
   mesh,
   dimensionedScalar("nEpsf",dimless, 0.0)
);

Info<< "create effective diffusivity field" << endl;
volScalarField Eps2Df
(
    "Eps2Df",
    eps*eps*Df
);
	//if (VoS=="VoS-psi")
    //{
		Info<<"creating effective specific surface area"<<endl;
		volScalarField ssArea
		(
			"ssArea",
			mag(fvc::grad(eps))*(4.0*eps*(1-eps))
		);
	//}
// rK.write();
// Foam::volScalarField K
// (
//     "K",
//     1./rK
// );
// K.write();
//#include "compressibleCreatePhi.H"

Info<< "calculating face flux field phi in porous medium" << endl;

surfaceScalarField phiByEpsf
(
    "phiByEpsf",
    phi*rEpsf
);
//volScalarField ByEps("ByEps",rEps);
Info<< "==Complete creating porous Fields== "<<nl<<endl;
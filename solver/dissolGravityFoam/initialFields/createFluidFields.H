
volScalarField HC
(
    IOobject
    (
        "HC",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh.C().component(2)
);	
volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	1 + beta0*(Y-Y0) 
);
// rhok.write();
volScalarField rhof
(
    IOobject
    (
        "rhof",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	rho0*rhok
);

Info<< "Reading field U" << endl;
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

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


Info<< "Reading field p_rgh" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
 //#include "createPhi.H"
  Info<< "Reading/calculating face flux field phi\n" << endl;
	 surfaceScalarField phi
	 (
		 IOobject
		 (
			 "phi",
			 runTime.timeName(),
			 mesh,
			 IOobject::READ_IF_PRESENT,
			 IOobject::NO_WRITE
		 ),
		 fvc::flux(U)
	 );
	 
  volScalarField p
 (
     IOobject
     (
         "p",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
     ),
     p_rgh + rhok*gh
 );
// Force p_rgh to be consistent with p
//p_rgh = p - rhof*gh;

label pRefCell = 0;
scalar pRefValue = 0.0;

 setRefCell
 (
     p,
     p_rgh,
     pimple.dict(),
     pRefCell,
     pRefValue
 );
 if (p_rgh.needReference())
 {
     p += dimensionedScalar
     (
         "p",
         p.dimensions(),
         pRefValue - getRefCellValue(p, pRefCell)
     );
 }
mesh.setFluxRequired(p_rgh.name());

dimensionedScalar deltaN
(
    "deltaN",
    dimensionSet(0, -1, 0, 0, 0, 0, 0),
    //solidProperties.lookup("deltaN")
	1e-13
);
//dimensionedScalar deltaN = 1e-13/pow(average(mesh.V()), 1.0/3.0);

Info<< "==Complete creating Fluid Fields== "<<nl<<endl;

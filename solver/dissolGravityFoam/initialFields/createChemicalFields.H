#include "dimensionSets.H"
#include "dimensionedScalarFwd.H"
Info<<"==Creating the chemical fields=="<<nl<<endl;
IOdictionary chemicalProperties
(
    IOobject
    (
        "chemicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar beta(chemicalProperties.lookup("beta"));
dimensionedScalar reactionRate(chemicalProperties.lookup("reactionRate"));
dimensionedScalar Yeq(chemicalProperties.lookup("Yeq"));
dimensionedScalar Y0(chemicalProperties.lookup("Y0"));
//dimensionedScalar reactionHeat(chemicalProperties.lookup("reactionHeat"));

Info<<"Stoichiometric Coefficient: "<<beta<<endl;
Info<<"const solid dissolution rate: "<<reactionRate<<endl;
Info<<"equilibrium mass concentration: "<<Yeq<<endl;
Info<<"initial mass concentration: "<<Y0<<endl;
//Info<<"reaction heat: "<<reactionHeat<<endl;

Info<<"reading reactant composition Y "<<endl;
volScalarField Y
(
    IOobject
    (
        "Y",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
/*volScalarField Yun
(
    IOobject
    (
        "Yun",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Yeq - Y
);*/

//sold dissoultion rate
volScalarField mdotSolid
(
    IOobject
    (
        "mdotSolid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Zero",dimensionSet(1,-3,-1,0,0,0,0),0.0)
);

//fluid species comsumption rate 
volScalarField mdotY
(
    IOobject
    (
        "mdotY",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Zero",dimensionSet(0,0,-1,0,0,0,0),0.0)
);




Info<< "==Complete Creating the chemical fields=="<<nl<<endl;

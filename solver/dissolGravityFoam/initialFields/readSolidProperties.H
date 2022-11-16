Info << "==Reading const solid thermo-physical transport properties==\n" << endl;

IOdictionary solidProperties
(
    IOobject
    (
        "solid",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar rhos
(
    "rhos",
    dimDensity,
    solidProperties.lookup("rhos")
);
Info<<"solid density: "<<rhos<<endl;

Info<< "==Complete reading solid Properties== "<<nl<<endl;

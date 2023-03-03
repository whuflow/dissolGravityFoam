
volScalarField rAU("rAU", 1.0/UEqn.A());
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

surfaceScalarField phig(-rAUf*ghf*fvc::snGrad(rhok)*mesh.magSf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::flux(HbyA)
      +rAUf*fvc::ddtCorr(rEps, U, phiByEpsf)  
      + phig
    )
);

constrainPressure(p_rgh, U, phiHbyA, rAUf);

// multiple the transient term by eps to ensure the DBS continuity equation
/*
fvScalarMatrix p_rghDDtEqn
(
  -  fvc::ddt(epsSolid) 
  + fvc::div(phiHbyA)
  ==
    fvm::Su(mdotSolid*(1.0-rho0/rhos),p_rgh)
);*/



    mdotY=ssArea*reactionRate;

	mdotSolid = mdotY*beta*(Yeq-Y);
	rhof = rho0*rhok;
	while (pimple.correctNonOrthogonal())
	{
		fvScalarMatrix p_rghEqn
		(
			  fvm::laplacian(rAUf, p_rgh)
		+ fvm::Su(mdotSolid*(1.0/rhof-1.0/rhos),p_rgh)
		==
		  fvc::div(phiHbyA)
	//	- fvc::ddt(epsSolid)
		);

    p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

    p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));
    //p_rghEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        // Calculate the conservative fluxes
        phi = phiHbyA -  p_rghEqn.flux(); 

        // Explicitly relax pressure for momentum corrector
        p_rgh.relax();

        // Correct the momentum source with the pressure gradient flux
        // calculated from the relaxed pressure
        U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
    }
}
    #include "continuityErrs.H"

    p = p_rgh + rhok*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rhok*gh;
    }

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Global
    solidEqn

Description
    Solve the continuity for density.

\*---------------------------------------------------------------------------*/


{
    Info << "solving solid equation"<<endl;
    // Foam::volScalarField effReactionRate("mdots")

	if (VoS=="VoS-psi")
    {
        mdotY=ssArea*reactionRate;
    }
    if (VoS=="iVoS")
    {
		volVectorField gradEps = fvc::grad(eps);
		volVectorField nEpsv = -gradEps/(mag(gradEps) + deltaN);//deltaN???
        volVectorField phiR = reactionRate*nEpsv;//nEpsf???

        volVectorField phiREps = phiR*eps;
        mdotY = eps*fvc::div(phiR) - fvc::div(phiREps);
    }
	mdotSolid = mdotY*beta*(Yeq-Y);
    
    fvScalarMatrix sEqn
    (
          fvm::ddt(epsSolid) 
        ==
        - fvm::Su(mdotSolid/rhos,epsSolid)
        + fvOptions(epsSolid)
    );

    sEqn.relax();
    fvOptions.constrain(sEqn);

    sEqn.solve();

    fvOptions.correct(epsSolid);

    epsSolid=max(min(epsSolid,1.0),0.0);
    eps=1.0-epsSolid;
    //update the geometrical  property of the porous medium 
	if (VoS=="VoS-psi")
    {
		ssArea=mag(fvc::grad(eps))*(4.0*eps*(1-eps));
	}

}

// ************************************************************************* //

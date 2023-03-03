
//update the filed related to porous medium
tEps=max(eps,minEps);
rK= rK0*(1.0-eps)*(1.0-eps)/(tEps*tEps*tEps);

rEps=1.0/tEps;
rEpsf=fvc::interpolate(rEps);
phiByEpsf = phi*rEpsf;

// Solve the Momentum equation
fvVectorMatrix UEqn
(
      rEps*(fvm::ddt(U) + fvm::div(phiByEpsf, U))
    - rEps*fvm::laplacian(nu, U)
    ==
    - fvm::SuSp(nu*rK,U)
    + fvOptions(U) 
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
        ==
        fvc::reconstruct
        (
            (
                - ghf*fvc::snGrad(rhok)
                - fvc::snGrad(p_rgh)
            )*mesh.magSf()
        )
    );

    fvOptions.correct(U);
}

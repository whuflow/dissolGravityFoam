{
    Eps2Df = eps*eps*Df;


     mdotY=ssArea*reactionRate;

	//C=Y/Yeq;
	//mdotC=mdotY/Yeq;Yun = (Yeq-Y);
    fvScalarMatrix YEqn
    (
		  fvm::ddt(eps, Y)
        + fvm::div(phi,Y)
        - fvm::laplacian(Eps2Df,Y)
        ==
		- fvm::Sp(mdotY,Y)     //implicit scheme
        + fvm::Su(mdotY*Yeq,Y) //explicit scheme                  //???????????????????
        + fvOptions(Y)
    );

    YEqn.relax();

    fvOptions.constrain(YEqn);

    YEqn.solve();

    fvOptions.correct(Y);
    Y.max(0.0);
    Y.min(Yeq);	//Y=max(min(Y,Yeq),0.0);
	//Yun = Yeq-Y;
    rhok = 1 + beta0*(Y-Y0);
	
}
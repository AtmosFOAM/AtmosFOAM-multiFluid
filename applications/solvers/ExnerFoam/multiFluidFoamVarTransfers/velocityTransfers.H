Info << "Velocity transfers based on vertical velocity variance:" << endl;
// The transfer terms are only formulated for nParts == 2
if (nParts != 2)
{
    FatalErrorIn("velocityTransfer.H")
        << " heat transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}

if ((wVarTransfer) and (!(directVarianceTransfer)))
//if (0)
{
    const dimensionedScalar sigmaRhoMin("",sigmaRho[0].dimensions(),scalar(1e-16));
    
    
    for(label ip = 0; ip < partNames.size(); ip++)
    {
        w[ip] = u[ip].component(vector::Z);
        Info << "w" << ip << " pre transfer: [" << min(w[ip]).value() << ", " << max(w[ip]).value() << "]" << endl;
        Info << "wVar" << ip << " pre Var transfer: [" << min(wVar[ip]).value() << ", " << max(wVar[ip]).value() << "]" << endl;
    }
    
    volScalarField wVarTransferRate = 0*transferRate[1] + 1 / wVarTimescale;
    if (localThetaVarTransfer)
    {
        wVarTransferRate *= dt * transferRate[1];
    }
    //Limit heat transfer to prevent instability.
    wVarTransferRate = min(wVarTransferRate, 1/dt);
    

    wVar[0] = max(wVar[0], 0*wVar[0]);
    wVar[1] = max(wVar[1], 0*wVar[1]);
    
    sigmaRhoPreTransfer = sigmaRho;
    wPreTransfer = w;
    wVarPreTransfer = wVar;

    wTransferTerm[0][1] = wVarTransferRate * 0.5 *
    (
        1 + Foam::erf
        (
            (w[0] - w[1]) / 
            max
            (
                wMin, 
                Foam::sqrt( 2 * (wVar[0] + wVar[1]) ) 
            )
        )
    );
    wTransferTerm[1][0] = wTransferTerm[0][1];

    wVarTransferTerm[0][1] = wVarTransferRate / Foam::sqrt
    ( 
        2*M_PI * max( wVarMin, wVar[0] + wVar[1] ) 
    ) * 
    Foam::exp
    (
        - 0.5 * sqr(w[1] - w[0]) / 
        max(wVarMin, wVar[0] + wVar[1])  
    );
    wVarTransferTerm[1][0] = wVarTransferTerm[0][1];

    if (varMassTransfer)
    {
        //Use variances to inform mass transfer between fluids.
        /*sigmaRho[0] += dt*
        (
            sigmaRhoPreTransfer[1]*wTransferTerm[1][0]
          - sigmaRhoPreTransfer[0]*wTransferTerm[0][1]
        );
        sigmaRho[1] += dt*
        (
            sigmaRhoPreTransfer[0]*wTransferTerm[0][1]
          - sigmaRhoPreTransfer[1]*wTransferTerm[1][0]
        );
        
        sigmaRho.updateSum();
        #include "sigma.H"*/
        
        massTransfer[0][1] = max(sigmaRhoPreTransfer[0]*wTransferTerm[0][1], massTransferZero);
        massTransfer[1][0] = max(sigmaRhoPreTransfer[1]*wTransferTerm[1][0], massTransferZero);
        #include "applyMassTransfer.H"
        sigmaRho.updateSum();
        transfer.updateSum();
        
        thetaPreTransfer = theta;
        thetaVarPreTransfer = thetaVar;
        thetaVarPreTransfer[0] += sqr(thetaPreTransfer[0]);
        thetaVarPreTransfer[1] += sqr(thetaPreTransfer[1]);

        theta[0] = ((1+transfer[0])*theta[0] + transfer[1]*theta[1])
                   /(1+transfer.sum());
        theta[1] = (theta[1] + transfer[0]*theta[0])/(1+transfer[0]);


        thetaVar[0] = ((1+transfer[0])*thetaVarPreTransfer[0] + transfer[1]*thetaVarPreTransfer[1])
                   /(1+transfer.sum()) - sqr(theta[0]);
        thetaVar[1] = ((1+transfer[1])*thetaVarPreTransfer[1] + transfer[0]*thetaVarPreTransfer[0])
                   /(1+transfer.sum()) - sqr(theta[1]);
        
        thetaVar[0] = max(thetaVar[0], 0*thetaVar[0]);
        thetaVar[1] = max(thetaVar[1], 0*thetaVar[1]);
        
        theta.updateSum();
        thetaVar.updateSum();
    }
    else
    {
        //Factor which prevents change in mass fractions.
        wTransferTerm[0][1] *= sigmaRho[1] / 
            max(sigmaRho[0]+sigmaRho[1], minSigmaDiv*(rho[0]+rho[1]));
        wTransferTerm[1][0] *= sigmaRho[0] / 
            max(sigmaRho[0]+sigmaRho[1], minSigmaDiv*(rho[0]+rho[1]));
        
        //Factor which prevents change in mass fractions.
        wVarTransferTerm[0][1] *= sigmaRho[1] / 
            max(sigmaRho[0]+sigmaRho[1], minSigmaDiv*(rho[0]+rho[1]));
        wVarTransferTerm[1][0] *= sigmaRho[0] / 
            max(sigmaRho[0]+sigmaRho[1], minSigmaDiv*(rho[0]+rho[1]));
    }

    

    //w transfer based on mean w
    w[0] *= sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
    w[0] -= dt*wTransferTerm[0][1]*wPreTransfer[0] * 
        sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
    w[0] += dt*wTransferTerm[1][0]*wPreTransfer[1] *
        sigmaRhoPreTransfer[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
    w[0] -= dt*wVarTransferTerm[0][1]*wVarPreTransfer[0] * 
        sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
    w[0] -= dt*wVarTransferTerm[1][0]*wVarPreTransfer[1] *
        sigmaRhoPreTransfer[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
        
    w[1] *= sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    w[1] += dt*wTransferTerm[0][1]*wPreTransfer[0] * 
        sigmaRhoPreTransfer[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    w[1] -= dt*wTransferTerm[1][0]*wPreTransfer[1] * 
        sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    w[1] += dt*wVarTransferTerm[0][1]*wVarPreTransfer[0] * 
        sigmaRhoPreTransfer[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    w[1] += dt*wVarTransferTerm[1][0]*wVarPreTransfer[1] *
        sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
        
    wVar[0] *= sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
    wVar[0] -= sqr(w[0]);
    wVar[0] += sqr(wPreTransfer[0]) *
        sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
    wVar[0] -= dt*wTransferTerm[0][1]*wVarPreTransfer[0] * 
        sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
    wVar[0] += dt*wTransferTerm[1][0]*wVarPreTransfer[1] *
        sigmaRhoPreTransfer[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
    wVar[0] -= dt*wTransferTerm[0][1]*sqr(wPreTransfer[0]) * 
        sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
    wVar[0] += dt*wTransferTerm[1][0]*sqr(wPreTransfer[1]) * 
        sigmaRhoPreTransfer[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
        
    wVar[0] -= 2*dt*wVarTransferTerm[0][1]*wVarPreTransfer[0] * wPreTransfer[0] *
        sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
    wVar[0] -= 2*dt*wVarTransferTerm[1][0]*wVarPreTransfer[1] * wPreTransfer[1] *
        sigmaRhoPreTransfer[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
    wVar[0] +=   dt*wVarTransferTerm[0][1]*wVarPreTransfer[0] * 
        (wPreTransfer[0] - wPreTransfer[1]) *
        wVarPreTransfer[0] / max(wVarPreTransfer[0] + wVarPreTransfer[1], wVarMin) *
        sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
    wVar[0] -=   dt*wVarTransferTerm[1][0]*wVarPreTransfer[1] * 
        (wPreTransfer[0] - wPreTransfer[1]) *
        wVarPreTransfer[1] / max(wVarPreTransfer[0] + wVarPreTransfer[1], wVarMin) *
        sigmaRhoPreTransfer[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);


    wVar[1] *= sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    wVar[1] -= sqr(w[1]);
    wVar[1] += sqr(wPreTransfer[1]) *
        sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    wVar[1] += dt*wTransferTerm[0][1]*wVarPreTransfer[0] * 
        sigmaRhoPreTransfer[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    wVar[1] -= dt*wTransferTerm[1][0]*wVarPreTransfer[1] * 
        sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    wVar[1] += dt*wTransferTerm[0][1]*sqr(wPreTransfer[0]) * 
        sigmaRhoPreTransfer[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    wVar[1] -= dt*wTransferTerm[1][0]*sqr(wPreTransfer[1]) * 
        sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
        
    wVar[1] += 2*dt*wVarTransferTerm[0][1]*wVarPreTransfer[0] * wPreTransfer[0] *
        sigmaRhoPreTransfer[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    wVar[1] += 2*dt*wVarTransferTerm[1][0]*wVarPreTransfer[1] * wPreTransfer[1] *
        sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    wVar[1] -=   dt*wVarTransferTerm[0][1]*wVarPreTransfer[0] * 
        (wPreTransfer[0] - wPreTransfer[1]) *
        wVarPreTransfer[0] / max(wVarPreTransfer[0] + wVarPreTransfer[1], wVarMin) * 
        sigmaRhoPreTransfer[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    wVar[1] +=   dt*wVarTransferTerm[1][0]*wVarPreTransfer[1] * 
        (wPreTransfer[0] - wPreTransfer[1]) *
        wVarPreTransfer[1] / max(wVarPreTransfer[0] + wVarPreTransfer[1], wVarMin) * 
        sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    

    /*volScalarField wTransferWvar1 = (w[1]-wPreTransfer[1]);
    volScalarField wTransferWvar0 = (w[0]-wPreTransfer[0]);
    volScalarField wVarTransferWvar1 = (wVar[1]-wVarPreTransfer[1]);
    volScalarField wVarTransferWvar0 = (wVar[0]-wVarPreTransfer[0]);
    
    w[0] = wPreTransfer[0] + wTransferWvar0 * pow( mag(wPreTransfer[1]/wOne-wPreTransfer[0]/wOne), 0.5);
    w[1] = wPreTransfer[1] + wTransferWvar1 * pow( mag(wPreTransfer[1]/wOne-wPreTransfer[0]/wOne), 0.5);
    wVar[0] = wVarPreTransfer[0] + wVarTransferWvar0 * pow( mag(wPreTransfer[1]/wOne-wPreTransfer[0]/wOne), 0.5);
    wVar[1] = wVarPreTransfer[1] + wVarTransferWvar1 * pow( mag(wPreTransfer[1]/wOne-wPreTransfer[0]/wOne), 0.5);*/
    
    // Update the sum
    w.updateSum();

    //Prevent rounding errors coursing negative variance.
    wVar[0] = max(wVar[0], 0*wVar[0]);
    wVar[1] = max(wVar[1], 0*wVar[1]);
    
    for(label ip = 0; ip < partNames.size(); ip++)
    {
        volVectorField du = (w[ip] - wPreTransfer[ip]) * Foam::vector(0,0,1);
        volFlux[ip] += (fvc::interpolate( du, "Uf") & mesh.Sf());
        sigmaFlux[ip] = sigmaRhof[ip]*volFlux[ip];
        u[ip] = fvc::reconstruct(volFlux[ip]);
        Uf[ip] = fvc::interpolate(u[ip], "Uf");
    }
    
    sigmaFlux.updateSum();
    volFlux.updateSum();
    u.updateSum();
    Uf.updateSum();

    for(label ip = 0; ip < nParts; ip++)
    {
        Info << "w" << ip << " post Var transfer: [" << min(w[ip]).value() << ", " << max(w[ip]).value() << "]" << endl;
    }
    for(label ip = 0; ip < nParts; ip++)
    {
        Info << "wVar" << ip << " post Var transfer: [" << min(wVar[ip]).value() << ", " << max(wVar[ip]).value() << "]" << endl;
    }
    Info << "TOTAL INTERNAL ENERGY POST TRANSFER: " << sum(w.sum()).value() << endl;
    
}



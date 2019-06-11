// advective form theta equation

for(label ip = 0; ip < nParts; ip++)
{
    thetaVar.ddt()[ip] = -fvc::div(volFlux[ip], thetaVar[ip], "theta")
                       + thetaVar[ip]*fvc::div(volFlux[ip]);
    
    // Viscous term
    if (alpha.value() > VSMALL)
    {
        thetaVar.ddt()[ip] += alpha * fvc::laplacian(thetaVar[ip]);
    }
    //thetaVar.ddt()[ip] += Ktheta * fvc::laplacian(thetaVar[ip]);
    
    if (wVarDiffusion)
    {
        thetaVar.ddt()[ip] += fvc::laplacian(KwVar[ip], thetaVar[ip]);
    }
    
    // Update theta according to the rate of change
    thetaVar[ip] = thetaVar[ip].oldTime() + dt*
    (
        (1-offCentre)*thetaVar.ddt()[ip].oldTime()
      + offCentre*thetaVar.ddt()[ip]
    );
    

    thetaVar[ip].correctBoundaryConditions();
}

// Mass transfer terms (for nParts = 2 only) block implicit solution
if (nParts != 2)
{
    FatalErrorIn("thetaEqn.H")
        << " transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}

thetaVar[0]  = max(thetaVar[0], 0*thetaVar[0]);
thetaVar[1]  = max(thetaVar[1], 0*thetaVar[1]);
thetaVarPreTransfer = thetaVar;

/*thetaVar[0] *= sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] -= sqr(theta[0]);
thetaVar[0] += sqr(thetaPreTransfer[0]) *
    sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] -= dt*transferRate[0]*thetaVarPreTransfer[0] * 
    sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] += dt*transferRate[1]*thetaVarPreTransfer[1] *
    sigmaRhoPreTransfer[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] -= dt*transferRate[0]*sqr(thetaPreTransfer[0]) * 
    sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] += dt*transferRate[1]*sqr(thetaPreTransfer[1]) * 
    sigmaRhoPreTransfer[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
    
thetaVar[1] *= sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] -= sqr(theta[1]);
thetaVar[1] += sqr(thetaPreTransfer[1]) *
    sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] += dt*transferRate[0]*thetaVarPreTransfer[0] * 
    sigmaRhoPreTransfer[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] -= dt*transferRate[1]*thetaVarPreTransfer[1] * 
    sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] += dt*transferRate[0]*sqr(thetaPreTransfer[0]) * 
    sigmaRhoPreTransfer[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] -= dt*transferRate[1]*sqr(thetaPreTransfer[1]) * 
    sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
*/

/*thetaVar[0] = ((1+transfer[0])*thetaVar[0] + transfer[1]*thetaVar[1])
           /(1+transfer.sum());
thetaVar[1] = (thetaVar[1] + transfer[0]*thetaVar[0])/(1+transfer[0]);
*/

/*thetaVarPreTransfer = thetaVar;
thetaVarPreTransfer[0] += sqr(thetaPreTransfer[0]);
thetaVarPreTransfer[1] += sqr(thetaPreTransfer[1]);

thetaVar[0] = ((1+transfer[0])*thetaVarPreTransfer[0] + transfer[1]*thetaVarPreTransfer[1])
           /(1+transfer.sum()) - sqr(theta[0]);
thetaVar[1] = ((1+transfer[1])*thetaVarPreTransfer[1] + transfer[0]*thetaVarPreTransfer[0])
           /(1+transfer.sum()) - sqr(theta[1]);*/


thetaVar[0] += sqr(thetaPreTransfer[0]); 
thetaVar[0] *= max(sigmaRhoPreTransfer[0], minSigmaDiv*rho[0]) / max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] -= dt*transferRate[0]*
    (thetaVarPreTransfer[0]+sqr(thetaPreTransfer[0])) * sigmaRhoPreTransfer[0] / 
    max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] += dt*transferRate[1]*
    (thetaVarPreTransfer[1]+sqr(thetaPreTransfer[1])) * sigmaRhoPreTransfer[1] / 
    max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] -= sigmaRhoPreTransfer[0] * 0.5 * 
    (thetaPreTransfer[0] + thetaTransferBound[0]) * 
    sqrt(0.5*thetaVarPreTransfer[0]/M_PI) * expTheta[0] / 
    max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] -= sigmaRhoPreTransfer[1] * 0.5 * 
    (thetaPreTransfer[1] + thetaTransferBound[1]) * 
    sqrt(0.5*thetaVarPreTransfer[1]/M_PI) * expTheta[1] / 
    max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] -= sqr(theta[0]);
thetaVar[0]  = max(thetaVar[0], 0*thetaVar[0]);

thetaVar[1] += sqr(thetaPreTransfer[1]);
thetaVar[1] *= max(sigmaRhoPreTransfer[1], minSigmaDiv*rho[1]) / max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] += dt*transferRate[0]*
    (thetaVarPreTransfer[0]+sqr(thetaPreTransfer[0])) * sigmaRhoPreTransfer[0] /
    max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] -= dt*transferRate[1]*
    (thetaVarPreTransfer[1]+sqr(thetaPreTransfer[1])) * sigmaRhoPreTransfer[1] /
    max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] += sigmaRhoPreTransfer[0] * 0.5 * 
    (thetaPreTransfer[0] + thetaTransferBound[0]) * 
    sqrt(0.5*thetaVarPreTransfer[0]/M_PI) * expTheta[0] /
    max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] += sigmaRhoPreTransfer[1] * 0.5 * 
    (thetaPreTransfer[1] + thetaTransferBound[1]) * 
    sqrt(0.5*thetaVarPreTransfer[1]/M_PI) * expTheta[1] /
    max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] -= sqr(theta[1]);
thetaVar[1]  = max(thetaVar[1], 0*thetaVar[1]);


/*
//thetaVar[0] *= 0;
//thetaVar[0] *= sigmaRhoMin/max(sigmaRho[0], sigmaRhoMin);
thetaVar[0] += sqr(thetaPreTransfer[0]); 
thetaVar[0] *= max(sigmaRhoPreTransfer[0], sigmaRhoMin) / max(sigmaRho[0], sigmaRhoMin);
thetaVar[0] -= dt*transferRate[0]*
    (thetaVarPreTransfer[0]+sqr(thetaPreTransfer[0])) * sigmaRhoPreTransfer[0] / 
    max(sigmaRho[0], sigmaRhoMin);
thetaVar[0] += dt*transferRate[1]*
    (thetaVarPreTransfer[1]+sqr(thetaPreTransfer[1])) * sigmaRhoPreTransfer[1] / 
    max(sigmaRho[0], sigmaRhoMin);
thetaVar[0] -= sigmaRhoPreTransfer[0] * 0.5 * 
    (thetaPreTransfer[0] + thetaTransferBound[0]) * 
    sqrt(0.5*thetaVarPreTransfer[0]/M_PI) * expTheta[0] / 
    max(sigmaRho[0], sigmaRhoMin);
thetaVar[0] -= sigmaRhoPreTransfer[1] * 0.5 * 
    (thetaPreTransfer[1] + thetaTransferBound[1]) * 
    sqrt(0.5*thetaVarPreTransfer[1]/M_PI) * expTheta[1] / 
    max(sigmaRho[0], sigmaRhoMin);
thetaVar[0] -= sqr(theta[0]);
thetaVar[0]  = max(thetaVar[0], 0*thetaVar[0]);

//thetaVar[1] *= 0;
//thetaVar[1] *= sigmaRhoMin/max(sigmaRho[1], sigmaRhoMin); 
thetaVar[1] += sqr(thetaPreTransfer[1]);
thetaVar[1] *= max(sigmaRhoPreTransfer[1], sigmaRhoMin) / max(sigmaRho[1], sigmaRhoMin);
thetaVar[1] += dt*transferRate[0]*
    (thetaVarPreTransfer[0]+sqr(thetaPreTransfer[0])) * sigmaRhoPreTransfer[0] /
    max(sigmaRho[1], sigmaRhoMin);
thetaVar[1] -= dt*transferRate[1]*
    (thetaVarPreTransfer[1]+sqr(thetaPreTransfer[1])) * sigmaRhoPreTransfer[1] /
    max(sigmaRho[1], sigmaRhoMin);
thetaVar[1] += sigmaRhoPreTransfer[0] * 0.5 * 
    (thetaPreTransfer[0] + thetaTransferBound[0]) * 
    sqrt(0.5*thetaVarPreTransfer[0]/M_PI) * expTheta[0] /
    max(sigmaRho[1], sigmaRhoMin);
thetaVar[1] += sigmaRhoPreTransfer[1] * 0.5 * 
    (thetaPreTransfer[1] + thetaTransferBound[1]) * 
    sqrt(0.5*thetaVarPreTransfer[1]/M_PI) * expTheta[1] /
    max(sigmaRho[1], sigmaRhoMin);
thetaVar[1] -= sqr(theta[1]);
thetaVar[1]  = max(thetaVar[1], 0*thetaVar[1]);*/

for(label ip = 0; ip < nParts; ip++)
{
    Info << "theta" << ip << " post Var transfer: [" << min(theta[ip]).value() << ", " << max(theta[ip]).value() << "]" << endl;
}
for(label ip = 0; ip < nParts; ip++)
{
    Info << "thetaVar" << ip << " post Var transfer: [" << min(thetaVar[ip]).value() << ", " << max(thetaVar[ip]).value() << "]" << endl;
}
Info << "TOTAL CHANGE IN INTERNAL ENERGY: " << sum(mag(theta.mean()-thetaPreTransfer.mean())).value() << endl;

thetaVar.updateSum();



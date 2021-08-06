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

/*thetaVarPreTransfer = thetaVar;

thetaVar[0] *= sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
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
thetaVar.updateSum();



// advective form w equation

for(label ip = 0; ip < nParts; ip++)
{
    wVar.ddt()[ip] = -fvc::div(volFlux[ip], wVar[ip], "theta")
                       + wVar[ip]*fvc::div(volFlux[ip]);
    
    // Viscous term
    if (nu.value() > VSMALL)
    {
        wVar.ddt()[ip] += nu * fvc::laplacian(wVar[ip]);
    }
    
    // Update theta according to the rate of change
    wVar[ip] = wVar[ip].oldTime() + dt*
    (
        (1-offCentre)*wVar.ddt()[ip].oldTime()
      + offCentre*wVar.ddt()[ip]
    );

    wVar[ip].correctBoundaryConditions();
}

// Mass transfer terms (for nParts = 2 only) block implicit solution
if (nParts != 2)
{
    FatalErrorIn("thetaEqn.H")
        << " transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}

wVarPreTransfer = wVar;

wVar[0] *= sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] -= sqr(w[0]);
wVar[0] += sqr(wPreTransfer[0]) *
    sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] -= dt*transferRate[0]*wVarPreTransfer[0] * 
    sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] += dt*transferRate[1]*wVarPreTransfer[1] *
    sigmaRhoPreTransfer[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] -= dt*transferRate[0]*sqr(wPreTransfer[0]) * 
    sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] += dt*transferRate[1]*sqr(wPreTransfer[1]) * 
    sigmaRhoPreTransfer[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
    
wVar[1] *= sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] -= sqr(w[1]);
wVar[1] += sqr(wPreTransfer[1]) *
    sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] += dt*transferRate[0]*wVarPreTransfer[0] * 
    sigmaRhoPreTransfer[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] -= dt*transferRate[1]*wVarPreTransfer[1] * 
    sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] += dt*transferRate[0]*sqr(wPreTransfer[0]) * 
    sigmaRhoPreTransfer[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] -= dt*transferRate[1]*sqr(wPreTransfer[1]) * 
    sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    
wVar.updateSum();



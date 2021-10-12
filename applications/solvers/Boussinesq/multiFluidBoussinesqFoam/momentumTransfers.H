// Velocity of the fluid transferred

for(label ip = 0; ip < nParts; ip++) for(label jp = ip+1; jp < nParts; jp++)
{
    volFluxT(ip,jp) = volFlux[ip];
    volFluxT(jp,ip) = volFlux[jp];
}

if (wTransferScale > SMALL && nParts == 2)
{
    volFluxT(0,1) += wTransferScale*sigmaf[0]*kdir*max(-volFlux[0]*kdir,fluxZero);
    volFluxT(1,0) -= wTransferScale*sigmaf[1]*kdir*max(volFlux[1]*kdir,fluxZero);
}
else if (wTransferScale > SMALL && nParts == 3)
{
    volFluxT(0,1) -= wTransferScale*sigmaf[0]*kdir*min(volFlux[0]*kdir,fluxZero);
    volFluxT(0,2) -= wTransferScale*sigmaf[0]*kdir*max(volFlux[0]*kdir,fluxZero);
    volFluxT(1,0) -= wTransferScale*sigmaf[1]*kdir*max(volFlux[1]*kdir,fluxZero);
    volFluxT(1,2) -= wTransferScale*sigmaf[1]*kdir*max(volFlux[1]*kdir,fluxZero);
    volFluxT(2,0) -= wTransferScale*sigmaf[2]*kdir*min(volFlux[2]*kdir,fluxZero);
    volFluxT(2,1) -= wTransferScale*sigmaf[2]*kdir*min(volFlux[2]*kdir,fluxZero);
/*
    surfaceScalarField transferTmp = sigmaf[0]*volFlux.ddt()[0]
                      + massTransferf(1,0)*(volFluxT(1,0)-volFlux[0])
                      + massTransferf(2,0)*(volFluxT(2,0)-volFlux[0]);
    volFluxT(0,1) += max(transferTmp, fluxChangeZero)
                    /max(massTransferf(0,1), rateSMALL);
    volFluxT(0,2) += min(transferTmp, fluxChangeZero)
                     max(massTransferf(0,2), rateSMALL);
*/
}

// Transfer the velocity
sigmaf.updateSum();
volFlux.transferField(interpolate(massTransferf, massTransfer), volFluxT, dt);

// Calculate and apply the drag between all pairs of fluids
if (mag(dragCoeff) > SMALL)
{
    // Calculate
    for(label ip = 0; ip < nParts; ip++) for(label jp=ip+1; jp < nParts; jp++)
    {
        surfaceScalarField dragCommon = dragCoeff/cloudRadius
                               * mag(fvc::interpolate(u[ip] - u[jp]));
        // Limit to avoid time step restrictions
        dragCommon = min(dragCommon, 1/dt);
        
        dragCommon *= volFlux[ip] - volFlux[jp];
        
        drag(ip,jp) = -sigmaf[jp]*dragCommon;
        drag(jp,ip) = sigmaf[ip]*dragCommon;
    }
    // Apply the drag
    for(label ip = 0; ip < nParts; ip++) for(label jp=ip+1; jp < nParts; jp++)
    {
        volFlux[ip] += dt*drag(ip,jp);
        volFlux[jp] += dt*drag(jp,ip);
    }
}


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
else if (wTransferScale > SMALL && nParts == 3) // Transfer dependent on TKE
{
    // Entrainmnet
    volFluxT(0,1) += wTransferScale*sqrt(linearInterpolate(2*TKE[0]))
                    *mesh.Sf().component(2);
    volFluxT(0,2) -= wTransferScale*sqrt(linearInterpolate(2*TKE[0]))
                    *mesh.Sf().component(2);

    // Detrainmnet
    volFluxT(1,0) -= wTransferScale*sqrt(linearInterpolate(2*TKE[1]))
                    *mesh.Sf().component(2);
    volFluxT(2,0) += wTransferScale*sqrt(linearInterpolate(2*TKE[2]))
                    *mesh.Sf().component(2);

    // trainmnet
    volFluxT(1,2) -= wTransferScale*sqrt(linearInterpolate(2*TKE[1]))
                    *mesh.Sf().component(2);
    volFluxT(2,1) += wTransferScale*sqrt(linearInterpolate(2*TKE[2]))
                    *mesh.Sf().component(2);
}


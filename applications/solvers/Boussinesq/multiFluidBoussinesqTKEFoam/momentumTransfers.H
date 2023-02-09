// Velocity of the fluid transferred

for(label ip = 0; ip < nParts; ip++) for(label jp = ip+1; jp < nParts; jp++)
{
    volFluxT(ip,jp) = volFlux[ip];
    volFluxT(jp,ip) = volFlux[jp];
}

if (wTransferScale > SMALL && nParts == 2)
{
    volFluxT(0,1) += wTransferScale*sqrt(linearInterpolate(2*TKE[0]))
                    *mesh.Sf().component(2);
    volFluxT(1,0) -= wTransferScale*sqrt(linearInterpolate(2*TKE[1]))
                    *mesh.Sf().component(2);
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
//    volFluxT(1,0) -= wTransferScale*sigmaf[1]*kdir*max(volFlux[1]*kdir,fluxZero);
//    volFluxT(2,0) += wTransferScale*sigmaf[2]*kdir*max(-volFlux[2]*kdir,fluxZero);

    // trainmnet
    volFluxT(1,2) -= wTransferScale*sqrt(linearInterpolate(2*TKE[1]))
                    *mesh.Sf().component(2);
    volFluxT(2,1) += wTransferScale*sqrt(linearInterpolate(2*TKE[2]))
                    *mesh.Sf().component(2);
//    volFluxT(1,2) -= wTransferScale*sigmaf[1]*kdir*max(volFlux[1]*kdir,fluxZero);
//    volFluxT(2,1) += wTransferScale*sigmaf[2]*kdir*max(-volFlux[2]*kdir,fluxZero);
}

// Calculate and apply the drag between all pairs of fluids
if (mag(dragRadius).value() > SMALL)
{
    // Calculate
    for(label ip = 0; ip < nParts; ip++) for(label jp=ip+1; jp < nParts; jp++)
    {
        surfaceScalarField dragCommon = 1/max(linearInterpolate(turbLength[0]), dragRadius)
                               * mag(fvc::interpolate(u[ip] - u[jp]));
        // Limit to avoid time step restrictions
        dragCommon = min(dragCommon, 1/dt);
        
        dragCommon *= volFlux[ip] - volFlux[jp];
        
        drag(ip,jp) = -sigmaf[jp]*dragCommon;
        drag(jp,ip) = sigmaf[ip]*dragCommon;
    }
    // Apply the drag and calculate the TKE sourse
    //dragTKEsource *= 0;
    for(label ip = 0; ip < nParts; ip++) for(label jp=ip+1; jp < nParts; jp++)
    {
        volFlux[ip] += dt*drag(ip,jp);
        volFlux[jp] += dt*drag(jp,ip);
        //dragTKEsource -= sigma[ip]*(u[ip] & fvc::reconstruct(drag(ip,jp)))
        //               + sigma[jp]*(u[jp] & fvc::reconstruct(drag(jp,ip)));
    }
}


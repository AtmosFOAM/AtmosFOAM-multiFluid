// Mass, buoyancy and momentum transfers due to diffusion

// Calculate and apply diffusive transfers
if (Ksigma.value() > SMALL)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        for(label jp=ip+1; jp < nParts; jp++)
        {
            volScalarField lapI = 2*Ksigma*fvc::laplacian(sigma[ip]);
            volScalarField lapJ = 2*Ksigma*fvc::laplacian(sigma[jp]);
        
            massTransfer(ip,jp)
                 = max(-lapI, rateZero)
                  *max( lapJ, rateZero)
                  /max(lapJ - lapI, rateSMALL);
            massTransfer(jp,ip)
                 = max(-lapJ, rateZero)
                  *max( lapI, rateZero)
                  /max(lapI - lapJ, rateSMALL);
        }
    }

    sigma.transferMass(massTransfer, dt);
    interpolate(sigmaf, sigma);

    // Diffusive transfers of buoyancy
    if (mag(Ksigma.value() - alpha.value()) > SMALL)
    {
        b.transferField(massTransfer, dt*(Ksigma-alpha)/Ksigma);
    }

    // Diffusive transfers of momentum
    if (mag(Ksigma.value() - nu.value()) > SMALL)
    {
        volFlux.transferField
        (
            interpolate(massTransferf, massTransfer),
            dt*(Ksigma-nu)/Ksigma
        );
    }
}


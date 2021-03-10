// Mass, buoyancy and momentum transfers due to diffusion

// Calculate and apply diffusive transfers
if (Ksigma.value() > SMALL)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        for(label jp=ip+1; jp < nParts; jp++)
        {
            volScalarField massTransferTmp = Ksigma/nParts*fvc::laplacian
            (
                sigma[ip] -sigma[jp]
            );
            massTransfer(ip,jp) = -min(massTransferTmp, rateZero);
            massTransfer(jp,ip) = -min(-massTransferTmp, rateZero);
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
            dt*(nu - alpha)/nu
        );
    }
}


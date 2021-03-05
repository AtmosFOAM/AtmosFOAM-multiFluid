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
            massTransferD[ip][jp] = -max(massTransferTmp, rateZero);
            massTransferD[jp][ip] = -max(-massTransferTmp, rateZero);
        }
    }

    for(label ip = 0; ip < nParts; ip++)
    {
        for(label jp=ip+1; jp < nParts; jp++)
        {
            sigma[ip] += dt*(massTransferD[jp][ip] - massTransferD[ip][jp]);
            sigma[jp] += dt*(massTransferD[ip][jp] - massTransferD[jp][ip]);
        }
    }
}

// Calculate sigma on the face
for(label ip = 0; ip < nParts; ip++)
{
    sigmaf[ip] = fvc::interpolate(sigma[ip], "sigma");
}

// Diffusive transfers of buoyancy
if (mag(Ksigma.value() - alpha.value()) > SMALL)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        for(label jp = ip+1; jp < nParts; jp++)
        {
            volScalarField bTransfer = dt*(Ksigma-alpha)/Ksigma*(b[jp] - b[ip]);

            b[ip] += bTransfer*massTransferD[jp][ip]/max(sigma[ip], minSigmaDiv);
            b[jp] -= bTransfer*massTransferD[ip][jp]/max(sigma[jp], minSigmaDiv);
        }
    }
}

// Diffusive transfers of momentum
if (mag(Ksigma.value() - nu.value()) > SMALL)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        for(label jp = ip+1; jp < nParts; jp++)
        {
            surfaceScalarField uTransfer = dt*(Ksigma-nu)/Ksigma
                                            *(volFlux[jp] - volFlux[ip]);

            volFlux[ip] += uTransfer*fvc::interpolate(massTransferD[jp][ip])
                          /max(sigmaf[ip], minSigmaDiv);
            volFlux[jp] -= uTransfer*fvc::interpolate(massTransferD[ip][jp])
                          /max(sigmaf[jp], minSigmaDiv);
        }
    }
}


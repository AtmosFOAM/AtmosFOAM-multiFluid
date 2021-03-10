// update sigma according to the continuity equation (without transfers)
for(label ip = 0; ip < nParts; ip++)
{
    sigma.ddt()[ip] = -fvc::div
    (
        ((1-offCentre)*volFlux[ip].oldTime() + offCentre*volFlux[ip]),
        sigma[ip].oldTime(),
        "sigma"
    );

    // Update sigma based on the rate of change
    sigma[ip] = sigma[ip].oldTime() + dt*sigma.ddt()[ip];
}

// Ensure that sum sigma = 1
sigma.updateSum();
for(label ip = 0; ip < nParts; ip++)
{
    sigma[ip] /= sigma.sum();
}
sigma.updateSum();

// Calculate, limit and apply diffusive transfers
if (Ksigma.value() > SMALL)
{
    // Calculate diffusive tranfer terms
    for(label ip = 0; ip < nParts; ip++)
    {
        massTransferD[ip] = -min(Ksigma*fvc::laplacian(sigma[ip]), rateZero);

        // Do not transfer out more than there exists
        massTransferD[ip] = min
        (
            massTransferD[ip],
            max(sigma[ip]-minSigma, dimlessZero)/dt
        );
    }
    for(label ip = 0; ip < nParts; ip++)
    {
        sigma[ip] += dt*(massTransferD[1-ip] - massTransferD[ip]);
    }
}

// Calculate sigma on the face
for(label ip = 0; ip < nParts; ip++)
{
    sigmaf[ip] = fvc::interpolate(sigma[ip], "sigma");
}

// Ensure that sum sigma = 1
sigmaf.updateSum();
for(label ip = 0; ip < nParts; ip++)
{
    sigmaf[ip] /= sigmaf.sum();
}
sigmaf.updateSum();


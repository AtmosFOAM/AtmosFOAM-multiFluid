// update sigma according to the continuity equation without transfers
for(label ip = 0; ip < nParts; ip++)
{
    // Advection
    sigma.ddt()[ip] = -fvc::div
    (
        ((1-offCentre)*volFlux[ip].oldTime() + offCentre*volFlux[ip]),
        sigma[ip].oldTime(),
        "sigma"
    );

    // Update sigma based on the rate of change
    sigma[ip] = sigma[ip].oldTime() + dt*sigma.ddt()[ip];
    sigma[ip].correctBoundaryConditions();
}

// Ensure that sum sigma = 1
sigma.updateSum();
for(label ip = 0; ip < nParts; ip++)
{
    sigma[ip] /= sigma.sum();
}
sigma.updateSum();

// Calculate sigma on the face
interpolate(sigmaf, sigma);


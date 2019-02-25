// update sigma according to the continuity equation (without transfers)
for(label ip = 0; ip < nParts; ip++)
{
    sigma.ddt()[ip] = -fvc::div
    (
        (1-offCentre)*volFlux[ip].oldTime() + offCentre*volFlux[ip],
        sigma[ip].oldTime(),
        "sigma"
    );

    // Update sigma based on the rate of change
    sigma[ip] = sigma[ip].oldTime() + dt*sigma.ddt()[ip];
}
// Correct so that the sigmas sum to 1
sigma.updateSum();
for(label ip = 0; ip < nParts; ip++)
{
    sigma[ip] /= sigma.sum();
    sigma.ddt()[ip] = (sigma[ip]- sigma[ip].oldTime())/dt;
}
sigma.updateSum();

// Update sigmaf
for(label ip = 0; ip < nParts; ip++)
{
    sigmaf[ip] = fvc::interpolate(sigma[ip]);
}
sigmaf.updateSum();
// Correct so that the sigmas sum to 1
sigmaf.updateSum();
for(label ip = 0; ip < nParts; ip++)
{
    sigmaf[ip] /= sigmaf.sum();
}
sigmaf.updateSum();


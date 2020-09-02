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

    sigmaf[ip] = fvc::interpolate(sigma[ip], "sigma");
}

// Ensure that sum sigma = 1
sigma.updateSum();
sigmaf.updateSum();
for(label ip = 0; ip < nParts; ip++)
{
    sigma[ip] /= sigma.sum();
    sigmaf[ip] /= sigmaf.sum();
}
sigma.updateSum();
sigmaf.updateSum();

// Correct sigma.ddt
for(label ip = 0; ip < nParts; ip++)
{
    sigma.ddt()[ip] = (sigma[ip] - sigma[ip].oldTime())/dt;
}


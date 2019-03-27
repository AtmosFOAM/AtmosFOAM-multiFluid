// update sigma according to the continuity equation (without transfers)
for(label ip = 0; ip < nParts; ip++)
{

    sigmaf[ip] = fvc::interpolate(sigma[ip].oldTime());

    sigma.ddt()[ip] = -fvc::div
    (
        ((1-offCentre)*volFlux[ip].oldTime() + offCentre*volFlux[ip])//,
        *sigmaf[ip]
        //sigma[ip].oldTime(),
        //"sigma"
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

// Correct sigma.ddt
for(label ip = 0; ip < nParts; ip++)
{
    sigma.ddt()[ip] = (sigma[ip] - sigma[ip].oldTime())/dt;
}
sigma.updateSum();
sigmaf.updateSum();


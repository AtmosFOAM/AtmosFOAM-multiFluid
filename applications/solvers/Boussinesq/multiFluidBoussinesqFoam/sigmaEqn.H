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
}

// Ensure that ddt(sigma[0]) = - ddt(sigma[1])
sigma.ddt()[0] = 0.5*(sigma.ddt()[0] - sigma.ddt()[1]);
sigma.ddt()[1] = -sigma.ddt()[0];

for(label ip = 0; ip < nParts; ip++)
{
    // Update sigma based on the rate of change
    sigma[ip] = sigma[ip].oldTime() + dt*sigma.ddt()[ip];
}
sigma.updateSum();
sigmaf.updateSum();


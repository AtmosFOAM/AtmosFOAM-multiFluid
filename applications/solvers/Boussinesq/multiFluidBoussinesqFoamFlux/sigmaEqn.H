// update sigma according to the continuity equation without transfers
for(label ip = 0; ip < nParts; ip++)
{
    // Advection
    sigma.ddt()[ip] = -fvc::div(flux[ip]);

    // Update sigma based on the rate of change
    sigma[ip] = sigma[ip].oldTime() + dt*
    (
        (1-offCentre)*sigma.ddt()[ip].oldTime()
      + offCentre*sigma.ddt()[ip]
    );
    sigma[ip].correctBoundaryConditions();
}

// Ensure that sum sigma = 1
sigma.updateSum();
for(label ip = 0; ip < nParts; ip++)
{
    sigma[ip] /= sigma.sum();
}
sigma.updateSum();

// Constraint to keep sigma => 0 and to ensure not more than all transferred
// in one time step (dt sigma S <= sigma)
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp = 0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            // Do not transfer out more than there exists
            massTransfer[ip][jp] = min
            (
                massTransfer[ip][jp],
                max(sigma[ip]-minSigma, dimlessZero)/dt
            );

            // Increase values that are <minSigma (shouldn't be needed)
            massTransfer[ip][jp] = max
            (
                massTransfer[ip][jp],
                max(minSigma-sigma[jp], dimlessZero)/dt
            );
        }
    }
}

// Transfers
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            sigma[ip] += dt*(massTransfer[jp][ip] - massTransfer[ip][jp]);
        }
    }
}

// Update boundary conditions
for(label ip = 0; ip < nParts; ip++)
{
    sigma[ip].correctBoundaryConditions();
}

// Calculate sigma on the face
for(label ip = 0; ip < nParts; ip++)
{
    sigmaf[ip] = fvc::interpolate(sigma[ip], "sigma");
}

// Ensure that sum sigma = 1
sigmaf.updateSum();
/*
for(label ip = 0; ip < nParts; ip++)
{
    sigmaf[ip] /= sigmaf.sum();
}
sigmaf.updateSum();
*/
Info << "sigma[0] goes from " << min(sigma[0]).value() <<  " to "
    << max(sigma[0]).value() << endl;


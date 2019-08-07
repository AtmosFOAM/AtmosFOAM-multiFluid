// advective form theta equation

for(label ip = 0; ip < nParts; ip++)
{
    theta.ddt()[ip] = -fvc::div(volFlux[ip], theta[ip], "theta")
                    + theta[ip]*fvc::div(volFlux[ip]);
    
    // Viscous term
    if (alpha.value() > VSMALL)
    {
        theta.ddt()[ip] += alpha * fvc::laplacian(theta[ip]);
    }
    
    // Update theta according to the rate of change
    theta[ip] = theta[ip].oldTime() + dt*
    (
        (1-offCentre)*theta.ddt()[ip].oldTime()
      + offCentre*theta.ddt()[ip]
    );

    theta[ip].correctBoundaryConditions();
}

// Mass transfer terms (for nParts = 2 only) block implicit solution
if (nParts != 2)
{
    FatalErrorIn("thetaEqn.H")
        << " transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}


theta.updateSum();

for(label ip = 0; ip < nParts; ip++)
{
    thetaf[ip] = fvc::interpolate(theta[ip]);
}
thetaf.updateSum();

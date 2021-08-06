// advective form theta equation

for(label ip = 0; ip < nParts; ip++)
{
    theta.ddt()[ip] = -fvc::div(volFlux[ip], theta[ip], "theta")
                    + theta[ip]*fvc::div(volFlux[ip]);

    // Update theta according to the rate of change
    theta[ip] = theta[ip].oldTime() + dt*theta.ddt()[ip];

    theta[ip].correctBoundaryConditions();
    thetaf[ip] = fvc::interpolate(theta[ip], "theta");
}
theta.updateSum();
thetaf.updateSum();

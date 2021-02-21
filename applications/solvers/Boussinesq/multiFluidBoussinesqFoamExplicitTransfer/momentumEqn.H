// Explicit part of velocity for each partition (without mass tranfers)
for(label ip = 0; ip < nParts; ip++)
{
    // Rate of change of volFlux without diffusion or pressure gradient
    volFlux.ddt()[ip] = (mesh.Sf() &
    (
      -fvc::interpolate
       (
           fvc::div(volFlux[ip], u[ip], "u") - u[ip]*fvc::div(volFlux[ip]),
           "convection"
       )
    ));
}

// Update volFlux
for(label ip = 0; ip < nParts; ip++)
{
    volFlux[ip] = volFlux[ip].oldTime()
                + dt*(1-offCentre)*volFlux.ddt()[ip].oldTime()
                + dt*offCentre*volFlux.ddt()[ip];
}

volFlux.updateSum();


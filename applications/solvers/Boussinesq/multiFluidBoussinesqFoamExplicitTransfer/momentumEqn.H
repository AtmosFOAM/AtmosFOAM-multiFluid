// Explicit part of velocity for each partition (without mass tranfers)
for(label ip = 0; ip < nParts; ip++)
{
    // Compact Laplacian of the velocity
    laplacianU 
        = (fvc::interpolate(fvc::laplacian(u[ip])) & mesh.Sf())
        + mesh.magSf()*(fvc::snGrad
            (
                fvc::div(volFlux[ip] - (linearInterpolate(u[ip]) & mesh.Sf()))
            )
          );

    // Rate of change of volFlux without pressure gradient
    volFlux.ddt()[ip] = (mesh.Sf() &
    (
      -fvc::interpolate
       (
           fvc::div(volFlux[ip], u[ip], "u") - u[ip]*fvc::div(volFlux[ip]),
           "convection"
       )
    )) + nu*laplacianU;
}

// Update volFlux
for(label ip = 0; ip < nParts; ip++)
{
    volFlux[ip] = volFlux[ip].oldTime()
                + dt*(1-offCentre)*volFlux.ddt()[ip].oldTime()
                + dt*offCentre*volFlux.ddt()[ip]
                - dt*offCentre*fvc::snGrad(p[ip])*mesh.magSf();
}

volFlux.updateSum();


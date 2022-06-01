// Compact diffusion of the velocity, common to all fluids
surfaceScalarField laplacianU
    = (fvc::interpolate(fvc::laplacian(nuTurb[0],u[0])) & mesh.Sf())
    + mesh.magSf()*nuTurb[0]*fvc::snGrad
    (
        fvc::div(volFlux[0] - (linearInterpolate(u[0]) & mesh.Sf()))
    );
for(label ip = 1; ip < nParts; ip++)
{
    laplacianU += (fvc::interpolate(fvc::laplacian(nuTurb[ip], u[ip])) & mesh.Sf())
            + mesh.magSf()*nuTurb[ip]*fvc::snGrad
            (
                fvc::div(volFlux[ip] - (linearInterpolate(u[ip]) & mesh.Sf()))
            );
}

// Explicit part of velocity for each partition (without mass tranfers)
for(label ip = 0; ip < nParts; ip++)
{
    // Rate of change of volFlux without pressure gradient
    volFlux.ddt()[ip] = (mesh.Sf() &
    (
      -fvc::interpolate
       (
           fvc::div(volFlux[ip], u[ip], "u") - u[ip]*fvc::div(volFlux[ip]),
           "convection"
       )
    )) + laplacianU
    - fvc::snGrad(p[ip])*mesh.magSf();
}

// Update volFlux
for(label ip = 0; ip < nParts; ip++)
{
    volFlux[ip] = volFlux[ip].oldTime()
                + dt*(1-offCentre)*volFlux.ddt()[ip].oldTime()
                + dt*offCentre*volFlux.ddt()[ip];
}

volFlux.updateSum();


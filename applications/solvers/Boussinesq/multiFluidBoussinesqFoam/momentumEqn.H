// Diffusive updates to the momentum
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

    volFlux.ddt()[ip] = nu*laplacianU;
    
    if (mag(nu - Ksigma).value() > SMALL)
    {
//        volFlux.ddt()[ip] += (Ksigma-nu)/Ksigma*(volFlux[1-ip] - volFlux[ip])
//          * fvc::interpolate(massTransferD[1-ip]/max(sigma[ip], minSigmaDiv));
        volFlux.ddt()[ip] += (Ksigma-nu)/Ksigma*
        (
            mesh.Sf()
          & fvc::interpolate
            (
                (u[1-ip] - u[ip])
               *massTransferD[1-ip]/max(sigma[ip], minSigmaDiv)
            )
        );
    }
}
for(label ip = 0; ip < nParts; ip++)
{
    volFlux[ip] = volFlux[ip].oldTime() + dt*volFlux.ddt()[ip];
}

// Explicit part of velocity for each partition (without mass tranfers)
for(label ip = 0; ip < nParts; ip++)
{
    // Rate of change of flux without diffusion or pressure gradient
    volFlux.ddt()[ip] = (mesh.Sf() &
    (
      -fvc::interpolate
       (
           fvc::div(volFlux[ip], u[ip], "u")
         - u[ip]*fvc::div(volFlux[ip]),
           "convection"
       )
    ))
  - fvc::snGrad(p[ip])*mesh.magSf();
}

for(label ip = 0; ip < nParts; ip++)
{
    volFlux[ip] += dt*(1-offCentre)*volFlux.ddt()[ip].oldTime()
                + dt*offCentre*volFlux.ddt()[ip];
}

volFlux.updateSum();


// Explicit part of velocity for each partition (without mass tranfers)
for(label ip = 0; ip < nParts; ip++)
{
    // Rate of change of flux without diffusion or pressure gradient
    flux.ddt()[ip] = (mesh.Sf() &
    (
      -fvc::interpolate
       (
           fvc::div(flux[ip], u[ip], "u"),
           "convection"
       )
    ))
    - sigmaf[ip]*fvc::snGrad(p[ip])*mesh.magSf();
}

// Update flux
for(label ip = 0; ip < nParts; ip++)
{
    flux[ip] = flux[ip].oldTime()
                + dt*(1-offCentre)*flux.ddt()[ip].oldTime()
                + dt*offCentre*flux.ddt()[ip];

    // Transfers
    for(label jp = 0; jp < nParts; jp++) if (ip != jp)
    {
        flux[ip] += dt*
        (
            fvc::interpolate(massTransfer[jp][ip])*volFluxT[jp]
          - fvc::interpolate(massTransfer[ip][jp])*volFluxT[ip]
        );
    }
}

flux.updateSum();

partitionedVolVectorField u1
(
    "u1",
    partNames,
    fvc::reconstruct(flux[0]/max(sigmaf[0], minSigmaDiv)),
    sigma
);
u1[1] = fvc::reconstruct(flux[1]/max(sigmaf[1], minSigmaDiv));
u1.write();


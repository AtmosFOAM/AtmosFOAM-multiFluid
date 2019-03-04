for (int corr=0; corr<nCorr; corr++)
{
    // Explicit part of velocity for each partition (without mass tranfers)
    for(label ip = 0; ip < nParts; ip++)
    {
        // Rate of change of flux without pressure gradient
        volFlux.ddt()[ip] = (mesh.Sf() &
        (
          -fvc::interpolate
           (
               fvc::div(volFlux[ip], u[ip], "u")
             - u[ip]*fvc::div(volFlux[ip])
             - nu*fvc::laplacian(u[ip]),
               "convection"
           )
        ))
        + fvc::interpolate(b[ip])*mesh.Sf().component(2);

        volFlux[ip] = volFlux[ip].oldTime()
                    + dt*(1-offCentre)*volFlux.ddt()[ip].oldTime()
                    + dt*offCentre*volFlux.ddt()[ip];
    }
    volFlux.updateSum();

    // Pressure equation
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix PEqn
        (
            fvc::div(volFlux.sum())
          - offCentre*dt*fvm::laplacian(sigma.sum(), P)
        );
        PEqn.setReference(0,0);
        PEqn.solve(P.name());
    }
    
    // Back-substitutions

    for(label ip = 0; ip < nParts; ip++)
    {
        volFlux.ddt()[ip] -= fvc::snGrad(P)*mesh.magSf();
        volFlux[ip] -= offCentre*dt*fvc::snGrad(P)*mesh.magSf();
    }

    // Include the mass transfer and drag implicitly (for nParts=2)
    transferf[0] = fvc::interpolate(transfer[0]);
    transferf[1] = fvc::interpolate(transfer[1]);
    transferf.updateSum();
    
    volFlux[0] =
    (
        (1+transferf[0])*volFlux[0] + transferf[1]*volFlux[1]
    )/(1+transferf.sum());
    volFlux[1] = (volFlux[1] + transferf[0]*volFlux[0])
                /(1+transferf[0]);

    // Update velocities based on the flux
    for(label ip = 0; ip < nParts; ip++)
    {
        u[ip] = fvc::reconstruct(volFlux[ip]);
    }
    volFlux.updateSum();
}
u.updateSum();


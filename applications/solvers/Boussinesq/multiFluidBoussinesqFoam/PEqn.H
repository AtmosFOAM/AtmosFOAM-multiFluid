for (int corr=0; corr<nCorr; corr++)
{
    // Explicit part of velocity for each partition (without mass tranfers)
    for(label ip = 0; ip < nParts; ip++)
    {
        // Compact vector Laplacian of the velocity
        Uf[ip] = linearInterpolate(u[ip]);
        surfaceScalarField laplacianU
            = fvc::snGrad(fvc::div(volFlux[ip]))*mesh.magSf()
            - (linearInterpolate(fvc::curl(fvc::curl(Uf[ip]))) & mesh.Sf());
        
        // Rate of change of flux without pressure gradient
        volFlux.ddt()[ip] = (mesh.Sf() &
        (
          -fvc::interpolate
           (
               fvc::div(volFlux[ip], u[ip], "u")
             - u[ip]*fvc::div(volFlux[ip]),
               "convection"
           )
        ))
      + nu*laplacianU;
//      - fvc::snGrad(Pi[ip])*mesh.magSf();

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
          + offCentre*dt*fvc::div(bf.sum())
          - offCentre*dt*fvm::laplacian(P)
        );
        PEqn.setReference(0,0);
        PEqn.solve();
    }
    
    // Back-substitutions
    for(label ip = 0; ip < nParts; ip++)
    {
        volFlux.ddt()[ip] += bf[ip] - fvc::snGrad(P)*mesh.magSf();
        volFlux[ip] += offCentre*dt*(bf[ip] - fvc::snGrad(P)*mesh.magSf());
    }
    volFlux.updateSum();

    // Update velocities based on the flux
    for(label ip = 0; ip < nParts; ip++)
    {
        u[ip] = fvc::reconstruct(volFlux[ip]);
    }
}



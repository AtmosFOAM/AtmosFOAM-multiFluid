// Separate pressure equation for each fluid

if (gamma.value() > VSMALL)
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

        volFlux[ip] = volFlux[ip].oldTime()
                    + dt*(1-offCentre)*volFlux.ddt()[ip].oldTime()
                    + dt*offCentre*volFlux.ddt()[ip];
    }

    volFlux.updateSum();
    
    // Store lagged mean pressure
    POld = P;
    for(label ip = 0; ip < nParts; ip++)
    {
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix PEqn
            (
                fvm::Sp(1/gamma,Pi[ip])
              + fvc::div(volFlux[ip])
              + offCentre*dt*fvc::div(bf[ip])
              - offCentre*dt*fvm::laplacian(Pi[ip])
              - offCentre*dt*fvc::laplacian(POld)
            );
            PEqn.solve();
        }
        
        // Back substitutions
        volFlux[ip] -= offCentre*dt*
                       (
                           fvc::snGrad(Pi[ip])*mesh.magSf()
                         + fvc::snGrad(POld)*mesh.magSf()
                         - bf[ip]
                       );
        volFlux.ddt()[ip] -= fvc::snGrad(Pi[ip])*mesh.magSf()
                             + fvc::snGrad(POld)*mesh.magSf()
                             - bf[ip];
        u[ip] = fvc::reconstruct(volFlux[ip]);
    }
}
volFlux.updateSum();

// Pressure equation for shared pressure

for (int corr=0; corr<nCorr; corr++)
{
    // Pressure equation
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix PEqn
        (
            fvc::div(volFlux.sum())
          - offCentre*dt*fvm::laplacian(P)
          + offCentre*dt*fvc::laplacian(POld)
        );
        PEqn.setReference(0,0);
        PEqn.solve();
    }
    
    // Back-substitutions
    for(label ip = 0; ip < nParts; ip++)
    {
        volFlux.ddt()[ip] += fvc::snGrad(POld)*mesh.magSf() - fvc::snGrad(P)*mesh.magSf();
        volFlux[ip] += offCentre*dt*(fvc::snGrad(POld)*mesh.magSf() - fvc::snGrad(P)*mesh.magSf());
    }
    volFlux.updateSum();

    // Update velocities based on the flux
    for(label ip = 0; ip < nParts; ip++)
    {
        u[ip] = fvc::reconstruct(volFlux[ip]);
    }
}

// Pressure gradient in each fluid including drag, for post processing
for(label ip = 0; ip < nParts; ip++)
{
    dPdz[ip] = mesh.Sf().component(2)/mesh.magSf()*
    (
        fvc::snGrad(P+Pi[ip])
      + (1-sigmaf[ip])*dragCommon*(volFlux[ip] - volFlux[1-ip])/mesh.magSf()
    );
}

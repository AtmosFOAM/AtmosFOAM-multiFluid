for (int corr=0; corr<nCorr; corr++)
{
    // Explicit part of velocity for each partition (without mass tranfers)
    for(label ip = 0; ip < nParts; ip++)
    {
        // Compact vector Laplacian of the velocity
        surfaceVectorField Uf = linearInterpolate(u[ip]);
        surfaceScalarField laplacianU
            = fvc::snGrad(fvc::div(volFlux[ip]))*mesh.magSf()
            - (linearInterpolate(fvc::curl(fvc::curl(Uf))) & mesh.Sf());
        
        // Rate of change of flux without pressure gradient
        volFlux.ddt()[ip] = (mesh.Sf() &
        (
          -fvc::interpolate
           (
               fvc::div(volFlux[ip], u[ip], "u")
             - u[ip]*fvc::div(volFlux[ip]),
               "convection"
           )
        )) + nu*laplacianU;

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
    
    // Pressure equations for each fluid separately
    #include "PiEqn.H"

    if (!noTransfers)
    {
        // Include the mass transfer and drag implicitly (for nParts=2)
        transferf[0] = fvc::interpolate(transfer[0], "transfer");
        transferf[1] = fvc::interpolate(transfer[1], "transfer");
        transferf.updateSum();
    
        // Velocity difference of fluid transferred
        volFluxt[0] = -volFlux[0]*mesh.Sf().component(2)/mesh.magSf();
        volFluxt[1] = -volFlux[1]*mesh.Sf().component(2)/mesh.magSf();
        
        volFlux[0] =
        (
            (1+transferf[0])*
            (
                volFlux[0]
              - dt*fvc::interpolate
                (
                    massTransfer[0][1]/sigma[0],
                    "massTransfer"
                )*volFluxt[0]
              + transferf[1]*volFluxt[1]
            )
          + transferf[1]*
            (
                volFlux[1]
              + transferf[0]*volFluxt[0]
              - dt*fvc::interpolate
                (
                    massTransfer[1][0]/sigma[1],
                    "massTransfer"
                )*volFluxt[1]
            )
        )/(1+transferf.sum());
        volFlux[1] = 
        (
            volFlux[1]
          + transferf[0]*volFlux[0]
          + transferf[0]*volFluxt[0]
          - dt*fvc::interpolate
            (
                massTransfer[1][0]/sigma[1],
                "massTransfer"
            )*volFluxt[1]
        )
        /(1+transferf[0]);
    }
    
    // Update velocities based on the flux
    for(label ip = 0; ip < nParts; ip++)
    {
        u[ip] = fvc::reconstruct(volFlux[ip]);
    }
    volFlux.updateSum();
}
u.updateSum();



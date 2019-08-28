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



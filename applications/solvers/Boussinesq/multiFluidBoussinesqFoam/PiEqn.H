// Separate pressure equation for each fluid

if (gamma.value() > VSMALL)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix PEqn
            (
//                fvm::Sp((1-sigma[ip])/gamma,Pi[ip])
//                fvm::Sp(1/(gamma*sqr(sigma[ip])),Pi[ip])
                fvm::Sp(1/gamma,Pi[ip])
              + fvc::div(volFlux[ip])
              - offCentre*dt*fvm::laplacian(Pi[ip])
            );
            PEqn.setReference(0,0);
            PEqn.solve();
        }
        
        // Back substitutions
        volFlux[ip] -= offCentre*dt*fvc::snGrad(Pi[ip])*mesh.magSf();

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

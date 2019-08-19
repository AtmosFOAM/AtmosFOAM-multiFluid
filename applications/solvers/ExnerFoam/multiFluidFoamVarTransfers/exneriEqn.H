//Separate pressure equation for each fluid.

if (Kw.value() > VSMALL)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix ExnerEqn
            (
                fvm::Sp(1/Kw, Exneri[ip])
              + fvc::div(volFlux[ip])
              - offCentre*dt*Cp*fvm::laplacian(theta[ip],Exneri[ip])
            );
            ExnerEqn.setReference(0,0);
            ExnerEqn.solve();
        }
        
        // Back substitutions
        volFlux[ip] -= offCentre*dt*Cp*thetaf[ip]*fvc::snGrad(Exneri[ip])*mesh.magSf();

        u[ip] = fvc::reconstruct(volFlux[ip]);
    }
}


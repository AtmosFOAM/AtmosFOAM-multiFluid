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
        Uf[ip] = fvc::interpolate(u[ip], "Uf");
    }
}
volFlux.updateSum();

//Diagnostic pressure gradient term for each fluid including drag
cloudRadius = max(cloudRadiusMin, cloudRadiusMax*sigmaf[0]*sigmaf[1]);
dragCommon = dragCoeff*sigmaRhof.sum()/cloudRadius * mag(Uf[0] - Uf[1]);
for(label ip = 0; ip < nParts; ip++)
{
    dExnerdz[ip] = mesh.Sf().component(2)/mesh.magSf()*
    (
        Cp*thetaf[ip]*fvc::snGrad(Exner + Exneri[ip])
      + dragCommon*sigmaf[1-ip]/fvc::interpolate(rho[ip],"rho") * 
            (volFlux[ip]-volFlux[1-ip])/mesh.magSf()
    );
}

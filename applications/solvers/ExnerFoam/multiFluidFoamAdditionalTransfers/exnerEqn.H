for (int corr=0; corr<nCorr; corr++)
{
    gradPcoeffSum = Cp*fvc::interpolate(theta.sum(), "theta");
    Psi = (1-kappa)/kappa*rho.sum()/Exner;
    
    // Explicit part of velocity for each partition (without mass tranfers)
    for(label ip = 0; ip < nParts; ip++)
    {
        // Rate of change of flux without gravity
        volFlux.ddt()[ip] = (mesh.Sf() &
        (
          -fvc::interpolate
           (
               fvc::div(volFlux[ip], u[ip], "u")
             - u[ip] * fvc::div(volFlux[ip]),
               "convection"
           )
        ))
        + gSf - Cp*thetaf[ip]*fvc::snGrad(Exner)*mesh.magSf();
        
        // Viscous term
        if (nu.value() > VSMALL)
        {
            volFlux.ddt()[ip] += mesh.Sf() & fvc::interpolate
            (
                nu * fvc::laplacian(u[ip]), "convection"
            );
        }
        
        // Parameterisation of pressure difference (i.e. deltaP[ip] = Kw * divU)
        volFlux.ddt()[ip] += mesh.magSf() * fvc::snGrad(Kw * fvc::div(volFlux[ip]));
        
        // Coriolis term
        if (magSqr(Omega.value()) > VSMALL)
        {
            volFlux.ddt()[ip] -= mesh.Sf() & (2*(Omega ^ Uf[ip]));
        }
        
        // Update the volFlux from the rate of change
        volFlux[ip] = volFlux[ip].oldTime()
                    + dt*(1-offCentre)*volFlux.ddt()[ip].oldTime()
                    + dt*offCentre*volFlux.ddt()[ip];
        
        // Mass Flux
        sigmaFlux[ip] = sigmaRhof[ip]*volFlux[ip];
    }
    sigmaFlux.updateSum();

    // Helmholtz equation
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        ExnerPrime *= 0;
        fvScalarMatrix ExnerEqn
        (
            fvm::Sp(Psi/dt, ExnerPrime)
          + fvc::ddt(rho.sum())
          + (1-offCentre)*fvc::div(sigmaFlux.sum().oldTime())
          + offCentre*fvc::div(sigmaFlux.sum())
          - fvm::laplacian(dt*sqr(offCentre)*gradPcoeffSum, ExnerPrime)
        );
        
        ExnerEqn.solve(Exner.name());
        Info << "Time " << runTime.timeName() << " ExnerPrime convergence "
             << (fvc::domainIntegrate(mag(ExnerPrime))/sum(mesh.V())).value()
             << endl;
    }
    
    // Back-substitutions
    Exner += ExnerPrime;
    #include "sigma.H"

    for(label ip = 0; ip < nParts; ip++)
    {
        volFlux[ip] -= dt*offCentre*Cp*thetaf[ip]
                       *fvc::snGrad(ExnerPrime)*mesh.magSf();

        volFlux.ddt()[ip] -= Cp*thetaf[ip]
                            *fvc::snGrad(ExnerPrime)*mesh.magSf();
    }
    
    //Transfer velocity
    if (useVolFluxTransfer)
    {
        #include "applyVolFluxTransfer.H"
    }
    else if (useFluxTransfer)
    {
        #include "applyFluxTransfer.H"
    }

    // Update mass flux and velocities based on the flux
    for(label ip = 0; ip < nParts; ip++)
    {
        sigmaFlux[ip] = sigmaRhof[ip]*volFlux[ip];
        u[ip] = fvc::reconstruct(volFlux[ip]);
        Uf[ip] = fvc::interpolate(u[ip], "Uf");
    }
    sigmaFlux.updateSum();
    volFlux.updateSum();
}
u.updateSum();
Uf.updateSum();


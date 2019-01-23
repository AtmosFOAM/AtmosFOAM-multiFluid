{
    for(label ip = 0; ip < nParts; ip++)
    {
        Psi[ip] = pow(rho[ip],(2*kappa-1)/(kappa-1))
                 *pow(R/pRef*theta[ip], kappa/(kappa-1));
    }
    Psi.updateSum();
    
    gradPcoeffSum = Cp*thetaf.sum();

    for (int corr=0; corr<nCorr; corr++)
    {
        // Explicit part of velocity for each partition (without mass tranfers)
        for(label ip = 0; ip < nParts; ip++)
        {
            // Rate of change of flux without gravity
            volFlux.ddt()[ip] = mesh.Sf() &
            (
              -fvc::interpolate
               (
                   fvc::div(volFlux[ip], u[ip], "u")
                 - u[ip] * fvc::div(volFlux[ip]),
                   "convection"
               )
            );
            
            // Viscous term
            if (nu.value() > VSMALL)
            {
                volFlux.ddt()[ip] += mesh.Sf() & fvc::interpolate
                (
                    nu * fvc::laplacian(u[ip]), "convection"
                );
            }
            
            // Coriolis term
            if (magSqr(Omega.value()) > VSMALL)
            {
                volFlux.ddt()[ip] -= mesh.Sf() & (2*(Omega ^ Uf[ip]));
            }
            
            // Flux with gravity
            volFlux[ip] = volFlux[ip].oldTime() + dt*
            (
                (1-offCentre)*volFlux.ddt()[ip].oldTime()
              + offCentre*(volFlux.ddt()[ip] + gSf)
            );
            sigmaFlux[ip] = sigmaRhof[ip]*volFlux[ip];
        }
        sigmaFlux.updateSum();
        // Needed for Exner boundary conditions
        surfaceScalarField rhog("rhog", sigmaRhof.sum()*gSf);

        // Helmholtz equation
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix ExnerEqn
            (
                fvm::ddt(Psi.sum(), Exner)
              + (1-offCentre)*fvc::div(sigmaFlux.sum().oldTime())
              + offCentre*fvc::div(sigmaFlux.sum())
              - fvm::laplacian(sqr(offCentre)*dt*gradPcoeffSum, Exner)
            );

            if(ucorr==nOuterCorr-1 && nonOrth==nNonOrthCorr && corr==nCorr-1)
            {
                ExnerEqn.solve(Exner.name() + "Final");
            }
            else
            {
                ExnerEqn.solve(Exner.name());
            }

            // Back-substitutions
            if (nonOrth == nNonOrthCorr)
            {
                for(label ip = 0; ip < nParts; ip++)
                {
                    volFlux[ip] -= offCentre*dt*Cp*thetaf[ip]
                                   *fvc::snGrad(Exner)*mesh.magSf();

                    volFlux.ddt()[ip] += gSf - Cp*thetaf[ip]
                                        *fvc::snGrad(Exner)*mesh.magSf();
                }

                // Include the mass transfer and drag implicitly (for nParts=2)
                transferf[0] = fvc::interpolate(transfer[0], "massTransfer")
                             + dt*dragCommon
                             *fvc::interpolate(sigma[0]/rho[1], "massTransfer");
                transferf[1] = fvc::interpolate(transfer[1], "massTransfer")
                             + dt*dragCommon
                             *fvc::interpolate(sigma[1]/rho[0], "massTransfer");
                transferf.updateSum();
                
                volFlux[0] =
                (
                    (1+transferf[0])*volFlux[0] + transferf[1]*volFlux[1]
                )/(1+transferf.sum());
                volFlux[1] = (volFlux[1] + transferf[0]*volFlux[0])
                            /(1+transferf[0]);

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
        }
    }
}
u.updateSum();
Uf.updateSum();


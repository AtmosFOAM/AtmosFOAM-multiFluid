for (int corr=0; corr<nCorr; corr++)
{
    // Solve Pressure equation for each fluid and back-substitute
    for(label ip = 0; ip < nParts; ip++)
    {
        // Density in this fluid
        rho[ip] = pRef/(R*theta[ip])*pow(Exner[ip], (1-kappa)/kappa);
        surfaceScalarField rhof = fvc::interpolate(rho[ip], "sigmaRho");

        // Calculate the drag
        surfaceScalarField drag = sigmaf[1-ip]/rhof*dragCommon
                                  *((Uf[ip] - Uf[1-ip]) & mesh.Sf());

        // Rate of change of flux
        volFlux.ddt()[ip] = (mesh.Sf() &
        (
          -fvc::interpolate
           (
               fvc::div(volFlux[ip], u[ip], "u")
             - u[ip] * fvc::div(volFlux[ip]),
               "convection"
           )
        ))
        + gSf - Cp*thetaf[ip]*fvc::snGrad(Exner[ip])*mesh.magSf() - drag;
        
        // Update the volFlux from the rate of change
        volFlux[ip] = volFlux[ip].oldTime()
                    + dt*volFlux.ddt()[ip];

        // Helmholtz equation
        volScalarField Psi("Psi", (1-kappa)/kappa*rho[ip]/Exner[ip]);
        surfaceScalarField gradPcoeff = Cp*rhof*thetaf[ip];
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix ExnerEqn
            (
                fvm::Sp(Psi/dt, ExnerPrime)
              + fvc::ddt(rho[ip])
              + fvc::div(rhof*volFlux[ip])
              - fvm::laplacian(dt*gradPcoeff, ExnerPrime)
            );
        
            ExnerEqn.solve("Exner");
            Info << "Time " << runTime.timeName() << " fluid " << ip
                 << " ExnerPrime convergence "
                 << (fvc::domainIntegrate(mag(ExnerPrime))/sum(mesh.V())).value()
                 << endl;
        }
    
        // Back-substitutions
        Exner[ip] += ExnerPrime;

        volFlux[ip] -= dt*Cp*thetaf[ip]
                       *fvc::snGrad(ExnerPrime)*mesh.magSf();

        volFlux.ddt()[ip] -= Cp*thetaf[ip]
                            *fvc::snGrad(ExnerPrime)*mesh.magSf();

        // Update mass flux and velocities based on the flux
        u[ip] = fvc::reconstruct(volFlux[ip]);
        Uf[ip] = fvc::interpolate(u[ip], "Uf");
    }
}
rho.updateSum();
volFlux.updateSum();
u.updateSum();
Uf.updateSum();


// Pressure equation for shared pressure
partitionedSurfaceScalarField volFluxPreTransfer = volFlux;
if (sharedDiffusion)
{
    // Compact Laplacian of the velocity
    laplacianU = (fvc::interpolate(fvc::laplacian(u.updateSum())) & mesh.Sf())
               + mesh.magSf()*(fvc::snGrad
                (
                    fvc::div(volFlux.updateSum()
                  - (linearInterpolate(u.sum()) & mesh.Sf()))
                ));
}

// Explicit part of velocity for each partition (without mass tranfers)
for(label ip = 0; ip < nParts; ip++)
{
    if (!sharedDiffusion)
    {
        // Compact Laplacian of the velocity
        laplacianU 
            = (fvc::interpolate(fvc::laplacian(u[ip])) & mesh.Sf())
            + mesh.magSf()*(fvc::snGrad
                (
                    fvc::div(volFlux[ip] - (linearInterpolate(u[ip]) & mesh.Sf()))
                )
              );
    }
    
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
  + 0*nu*laplacianU
  //+ ( mesh.Sf() & fvc::interpolate((1/sigma[ip])*fvc::laplacian(sigma[ip]*nu, u[ip]), "convection" ) )
  - 0*fvc::snGrad(Pi[ip] - Pi.sum())*mesh.magSf();

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
forAll (volFlux[0], celli)
{
    if (mag(volFlux[0][celli]-volFluxPreTransfer[0][celli]) > 1e-8)
    {
        Info << endl;
        Info << "stage Peqn " << endl;
        Info << "t=" << runTime.timeName() << ", ucorr=" << ucorr << ", corr=" << corr << " " << endl;
        Info << celli << endl;
        Info << "sigma0: " << sigmaf[0][celli] << endl;
        Info << "sigma1: " << sigmaf[1][celli] << endl;
        Info << "S0: " << Sf[0][celli] << endl;
        Info << "S1: " << Sf[1][celli] << endl;
        Info << "w0: " << volFluxPreTransfer[0][celli]/2e7 << " to " << volFlux[0][celli]/2e7 << endl;
        Info << "w1: " << volFluxPreTransfer[1][celli]/2e7 << " to " << volFlux[1][celli]/2e7 << endl;
        Info << "w.ddt0: " << volFlux.ddt()[0][celli]/2e7 << endl;
        Info << "w.ddt1: " << volFlux.ddt()[1][celli]/2e7 << endl;
        
    }
}
// Back-substitutions
for(label ip = 0; ip < nParts; ip++)
{
    volFlux.ddt()[ip] += bf[ip] - fvc::snGrad(P)*mesh.magSf();
    volFlux[ip] += offCentre*dt*(bf[ip] - fvc::snGrad(P)*mesh.magSf());
}
volFlux.updateSum();



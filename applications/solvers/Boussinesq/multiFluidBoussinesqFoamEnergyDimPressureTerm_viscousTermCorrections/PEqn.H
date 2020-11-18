// Pressure equation for shared pressure

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
  + nu*laplacianU
  /*+ 2*nu*(fvc::interpolate
    (
        (fvc::grad(log(max(sigma[ip], VSMALL))) & fvc::grad(u[ip])),
        "convection"
    ) & mesh.Sf())
  + nu*volFlux[ip]*fvc::interpolate
    (
        (1/max(sigma[ip],VSMALL))*fvc::laplacian(sigma[ip]),
        "convection"
    )
  - 2*nu*(fvc::interpolate
    (
        (fvc::grad(log(max(sigma[ip], VSMALL))) & fvc::grad(u.sum())),
        "convection"
    ) & mesh.Sf())
  - nu*volFlux.sum()*fvc::interpolate
    (
        (1/max(sigma[ip],VSMALL))*fvc::laplacian(sigma[ip]),
        "convection"
    )*/
  - fvc::snGrad(Pi[ip])*mesh.magSf();// - Pi.sum())*mesh.magSf();
  //- fvc::interpolate(Pi[ip] - Pi.sum())*fvc::snGrad(log(max(sigma[ip], VSMALL)))*mesh.magSf(); // pressure gradient correction: allows sum of pressure terms = grad(P); allows for definite energy sink

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
volFlux.updateSum();

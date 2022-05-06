// flux form TKE equation without diffusion or transfers (explicit solution)
for(label ip = 0; ip < nParts; ip++)
{
    u[ip] = fvc::reconstruct(volFlux[ip]);
    volTensorField gradu = fvc::grad(u[ip]);
    volScalarField dbdz = fvc::div(bf[ip]*magk);
    volScalarField nuTurbC = turbLength[ip]*sqrt(TKE[ip]);
    
    TKE.ddt()[ip] = - fvc::div(volFlux[ip], TKE[ip], "TKE")
                  + TKE[ip]*fvc::div(volFlux[ip])
                  + fvc::laplacian(nuTurb[ip], TKE[ip])
                  + (nuTurbC + nu)*(dev(twoSymm(gradu)) && gradu)
                  - nuTurbC*dbdz
                  - 0.09*pow(TKE[ip], 1.5)/turbLength[ip];
    // Limit to maintain positivity
    TKE.ddt()[ip] = max(TKE.ddt()[ip], -(0.5-SMALL)*TKE[ip].oldTime()/dt);
    TKE.ddt()[ip].oldTime() 
        = max(TKE.ddt()[ip].oldTime(), -(0.5-SMALL)*TKE[ip].oldTime()/dt);

    // Update TKE according to rate of change
    TKE[ip] = TKE[ip].oldTime() + dt*
    (
        (1-offCentre)*TKE.ddt()[ip].oldTime()
      + offCentre*TKE.ddt()[ip]
    );

    TKE[ip].correctBoundaryConditions();
    
    // Update nuTurb
    nuTurb[ip] = nu + linearInterpolate(turbLength[ip]*sqrt(TKE[ip]));
}

TKE.updateSum();


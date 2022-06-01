// Diffusion common for all fluids
volScalarField kDiffusion = sigma[0]*fvc::laplacian(nuTurb[0], TKE[0]);
for(label ip = 1; ip < nParts; ip++)
{
    kDiffusion += sigma[ip]*fvc::laplacian(nuTurb[ip], TKE[ip]);
}

// TKE equation without diffusion or transfers (explicit solution)
for(label ip = 0; ip < nParts; ip++)
{
    u[ip] = fvc::reconstruct(volFlux[ip]);
    volTensorField gradu = fvc::grad(u[ip]);
    volScalarField dbdz = fvc::div(bf[ip]*magk);
    volScalarField nuTurbC = turbLength[ip]*sqrt(TKE[ip]) + nu;
    
    TKE.ddt()[ip] = - fvc::div(volFlux[ip], TKE[ip], "TKE")
                  + TKE[ip]*fvc::div(volFlux[ip])
                  + kDiffusion
                  + nuTurbC*(dev(twoSymm(gradu)) && gradu)
                  - nuTurbC*dbdz
                  - Cmu*pow(TKE[ip], 1.5)/turbLength[ip]
                  + gamma*(u[ip] & fvc::laplacian(u[ip]));

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
//    Info << "In kEqn.H TKE[" << ip << "] goes from " << min(TKE[ip]).value()
//         << " to " << max(TKE[ip]).value() << endl;
    
    // Update nuTurb
    nuTurb[ip] = nu + linearInterpolate(turbLength[ip]*sqrt(TKE[ip]));
}

TKE.updateSum();


// update rhoSigma according to the continuity equation (without transfers)
for(label ip = 0; ip < nParts; ip++)
{
    // Calculate the rate of change using old sigmaRho for bounded advection
    flux[ip] = linearInterpolate(rho[ip])*volFlux[ip];
    
    sigmaRho.ddt()[ip] = -fvc::div
    (
        (1-offCentre)*flux[ip].oldTime() + offCentre*flux[ip],
        sigma[ip].oldTime(),
        "sigmaRho"
    );
    
    //sigmaRho.ddt()[ip] += Ktheta * fvc::laplacian(sigmaRho[ip]);
    

    // Update sigmaRho based on the rate of change
    sigmaRho[ip] = sigmaRho[ip].oldTime() + dt*sigmaRho.ddt()[ip];
    
    
    
    sigmaRho[ip].correctBoundaryConditions();
}
flux.updateSum();
sigmaRho.updateSum();


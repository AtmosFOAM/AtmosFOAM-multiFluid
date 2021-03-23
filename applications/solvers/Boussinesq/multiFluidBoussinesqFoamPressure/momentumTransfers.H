// Velocity of the fluid transferred

for(label ip = 0; ip < nParts; ip++) for(label jp = ip+1; jp < nParts; jp++)
{
    volFluxT(ip,jp) = volFlux[ip];
    volFluxT(jp,ip) = volFlux[jp];
}

if (mag(wTransfer0) > SMALL)
{
    // Velocity anomalies of fluid transferred
    
    if (nParts >= 2)
    {
        volFluxT(0,1) += sigmaf[0]*wTransfer0
                       *mag(volFlux[0]*mesh.Sf().component(2))
                       *mesh.Sf().component(2)/sqr(mesh.magSf());
        volFluxT(1,0) -= sigmaf[1]*wTransfer0
                       *mag(volFlux[1]*mesh.Sf().component(2))
                       *mesh.Sf().component(2)/sqr(mesh.magSf());
    }
    if (nParts >= 3)
    {
        volFluxT(0,2) -= sigmaf[0]*wTransfer0
                       *mag(volFlux[0]*mesh.Sf().component(2))
                       *mesh.Sf().component(2)/sqr(mesh.magSf());
        volFluxT(2,0) += sigmaf[2]*wTransfer0
                       *mag(volFlux[2]*mesh.Sf().component(2))
                       *mesh.Sf().component(2)/sqr(mesh.magSf());
        volFluxT(1,2) -= sigmaf[1]*wTransfer0
                        *mag(volFlux[1]*mesh.Sf().component(2))
                        *mesh.Sf().component(2)/sqr(mesh.magSf());
        volFluxT(2,1) += sigmaf[2]*wTransfer0
                        *mag(volFlux[2]*mesh.Sf().component(2))
                        *mesh.Sf().component(2)/sqr(mesh.magSf());
    }
}

// Transfer the velocity
volFlux.transferField(interpolate(massTransferf, massTransfer), volFluxT, dt);

// Calculate and apply the drag between all pairs of fluids
if (mag(dragCoeff) > SMALL)
{
    // Calculate
    for(label ip = 0; ip < nParts; ip++) for(label jp=ip+1; jp < nParts; jp++)
    {
        surfaceScalarField dragCommon = dragCoeff/cloudRadius
                               * mag(fvc::interpolate(u[ip] - u[jp]));
        // Limit to avoid time step restrictions
        dragCommon = min(dragCommon, 1/dt);
        
        dragCommon *= volFlux[ip] - volFlux[jp];
        
        drag(ip,jp) = -sigmaf[jp]*dragCommon;
        drag(jp,ip) = sigmaf[ip]*dragCommon;
    }
    // Apply the drag
    for(label ip = 0; ip < nParts; ip++) for(label jp=ip+1; jp < nParts; jp++)
    {
        volFlux[ip] += dt*drag(ip,jp);
        volFlux[jp] += dt*drag(jp,ip);
    }
}

// Update velocities based on the volFlux
for(label ip = 0; ip < nParts; ip++)
{
    u[ip] = fvc::reconstruct(volFlux[ip]);
    divu[ip] = fvc::div(volFlux[ip]);
}
u.updateSum();
divu.updateSum();


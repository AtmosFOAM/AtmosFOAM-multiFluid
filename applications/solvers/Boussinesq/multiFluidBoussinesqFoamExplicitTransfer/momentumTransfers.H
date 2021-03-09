// Velocity of the fluid transferred

for(label ip = 0; ip < nParts; ip++) for(label jp = ip+1; jp < nParts; jp++)
{
    volFluxT(ip,jp) = volFlux[ip];
    volFluxT(jp,ip) = volFlux[jp];
}

if (wTransfer0)
{
    // only formulated for nParts == 2
    if (nParts != 2)
    {
        FatalErrorIn("momentumTransfer.H")
            << " wTransfer0 only formulated for nParts = 2, not nParts = "
            << nParts << exit(FatalError);
    }

    // Velocity anomalies of fluid transferred
    volFluxT(0,1) += sigmaf[0]*wPrimeCoeff
                              *mag(volFlux[0]*mesh.Sf().component(2))
                              *mesh.Sf().component(2)/sqr(mesh.magSf());
    volFluxT(1,0) -= sigmaf[1]*wPrimeCoeff
                              *mag(volFlux[1]*mesh.Sf().component(2))
                              *mesh.Sf().component(2)/sqr(mesh.magSf());
}

// Transfer the velocity
volFlux.transferField(interpolate(massTransferf, massTransfer), volFluxT, dt);

// Calculate and apply the drag
if (mag(dragCoeff) > SMALL)
{
    // This code only works for nParts == 2
    if (nParts != 2)
    {
        FatalErrorIn("momentumTransfer.H")
            << " drag is only formulated for nParts = 2, not nParts = "
            << nParts << exit(FatalError);
    }

    // Calculate the cloud radius and hence the drag
    surfaceScalarField drag = dt*dragCoeff/cloudRadius
        * mag(fvc::interpolate(u[0] - u[1]));
    // Limit and multiply by flux difference:
    drag = min(drag, scalar(1));
    drag *= volFlux[0] - volFlux[1];
    
    volFlux[0] -= sigmaf[1]*drag;
    volFlux[1] += sigmaf[0]*drag;
}

// Update velocities based on the volFlux
for(label ip = 0; ip < nParts; ip++)
{
    u[ip] = fvc::reconstruct(volFlux[ip]);
}
u.updateSum();


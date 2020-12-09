// Updates due to diffusion (1st order in time). Applied first so that b.ddt
// can be re-used
for(label ip = 0; ip < nParts; ip++)
{
    b.ddt()[ip] = alpha*fvc::laplacian(b[ip]);

    if (mag(alpha - Ksigma).value() > SMALL)
    {
        b.ddt()[ip] += (Ksigma-alpha)/Ksigma*(b[1-ip] - b[ip])
            *massTransferD[1-ip]/max(sigma[ip], minSigmaDiv);
    }
}
for(label ip = 0; ip < nParts; ip++)
{
    b[ip] = b[ip].oldTime() + dt*b.ddt()[ip];
    b[ip].correctBoundaryConditions();
}

// advective form b equation without diffusion
for(label ip = 0; ip < nParts; ip++)
{
    b.ddt()[ip] = -fvc::div(volFlux[ip], b[ip], "b")
                + b[ip]*fvc::div(volFlux[ip]);

    if (!adjustQ) b.ddt()[ip] += Q;
}

for(label ip = 0; ip < nParts; ip++)
{
    // Update b according to the rate of change
    b[ip] += dt*
    (
        (1-offCentre)*b.ddt()[ip].oldTime()
      + offCentre*b.ddt()[ip]
    );
    b[ip].correctBoundaryConditions();
}

if (adjustQ)
{
    dimensionedScalar totalHeat("totalHeat", b[0].dimensions(), scalar(0));
    for(label ip = 0; ip < nParts; ip++)
    {
        totalHeat += fvc::domainIntegrate(sigma[ip]*b[ip])/sum(mesh.V());
    }
    Q = (totalHeat0 - totalHeat)/dt;
    for(label ip = 0; ip < nParts; ip++)
    {
        b[ip] += dt*Q;
    }
}

// bouyancy force on faces
for(label ip = 0; ip < nParts; ip++)
{
    bf[ip] = fvc::interpolate(b[ip], "b")*mesh.Sf().component(2);
}

b.updateSum();
bf.updateSum();


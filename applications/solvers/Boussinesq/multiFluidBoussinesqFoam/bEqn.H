// advective form b equation

for(label ip = 0; ip < nParts; ip++)
{
    b.ddt()[ip] = -Nsquared*u[ip].component(2)
                - fvc::div(volFlux[ip], b[ip], "b")
                + b[ip]*fvc::div(volFlux[ip])
                + alpha*fvc::laplacian(b[ip]);
}

if (adjustQ)
{
    Q = dimensionedScalar("zeroQ", Q.dimensions(), scalar(0));
    for(label ip = 0; ip < nParts; ip++)
    {
        Q -= fvc::domainIntegrate(sigma[ip]*b.ddt()[ip])/sum(mesh.V());
    }
}

for(label ip = 0; ip < nParts; ip++)
{
    b.ddt()[ip] += Q;
}

for(label ip = 0; ip < nParts; ip++)
{
    // Update b according to the rate of change
    b[ip] = b[ip].oldTime() + dt*
    (
        (1-offCentre)*b.ddt()[ip].oldTime()
      + offCentre*b.ddt()[ip]
    );

    b[ip].correctBoundaryConditions();
}

// Mass transfer terms (for nParts = 2 only) block implicit solution
if (nParts != 2)
{
    FatalErrorIn("bEqn.H")
        << " transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}

b[0] = ((1+transfer[0])*b[0] + transfer[1]*b[1])
       /(1+transfer.sum());

b[1] = (b[1] + transfer[0]*b[0])
       /(1+transfer[0]);

// bouyancy force on faces
for(label ip = 0; ip < nParts; ip++)
{
    bf[ip] = fvc::interpolate(b[ip], "b")*mesh.Sf().component(2);
}

b.updateSum();
bf.updateSum();


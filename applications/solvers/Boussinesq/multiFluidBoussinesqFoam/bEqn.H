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

// Bouyancy anomalies of fluid transferred
if (transferAnomalies)
{
    bt[0] = min
    (
        -sigma[0]*bTransferTimeScale*bTransferDiffusivity*min
        (
            fvc::laplacian(b[0]),
            dimensionedScalar("", dimensionSet(0,-1,-2,0,0), scalar(0))
        ),
        bTransferMax
    );

    bt[1] = max
    (
        -sigma[1]*bTransferTimeScale*bTransferDiffusivity*max
        (
            fvc::laplacian(b[1]),
            dimensionedScalar("", dimensionSet(0,-1,-2,0,0), scalar(0))
        ),
        -bTransferMax
    );
}

b[0] = (
           (1+transfer[0])*
           (
               b[0]
             - dt*massTransfer[0][1]/sigma[0]*bt[0]
             + transfer[1]*bt[1]
           )
         + transfer[1]*
           (
               b[1]
             - dt*massTransfer[1][0]/sigma[1]*bt[1]
             + transfer[0]*bt[0]
           )
       )
       /(1+transfer.sum());

b[1] = (
           b[1]
         + transfer[0]*(b[0] + bt[0])
         - dt*massTransfer[1][0]/sigma[1]*bt[1]
       )
       /(1+transfer[0]);

// bouyancy force on faces
for(label ip = 0; ip < nParts; ip++)
{
    bf[ip] = fvc::interpolate(b[ip], "b")*mesh.Sf().component(2);
}

b.updateSum();
bf.updateSum();


// advective form b equation

for(label ip = 0; ip < nParts; ip++)
{
    b.ddt()[ip] = -fvc::div(volFlux[ip], b[ip], "b")
                + b[ip]*fvc::div(volFlux[ip])
                + alpha*fvc::laplacian(b[ip]);
}

if (!adjustQ)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        b.ddt()[ip] += Q;
    }
}

if (mag(alpha - Ksigma).value() > SMALL)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        b.ddt()[ip] += (alpha - Ksigma)/Ksigma*(b[ip] - b[1-ip])
            *transferD[1-ip];
    }
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

if (!noTransfers)
{
    // Bouyancy anomalies of fluid transferred
    // Bouyancy transferred is linear combination of buoyancy in each fluid
    if (bAnomalyLinearCoeff > SMALL)
    {
        bt[0] = bAnomalyLinearCoeff*(b[1] - b[0]);
        bt[1] = - bt[0];
    }
    
    // Bouyancy transferred is constant factor times mean buoyancy
    if (bAnomalyMeanCoeff > SMALL)
    {
        bt[0] = bAnomalyMeanCoeff*sigma[0]*mag(b[0]);
        bt[1] = -bAnomalyMeanCoeff*sigma[1]*mag(b[1]);
    }

    // Transfers
    b[0] = (
               (1+transfer[0])*(b[0] - S[0]*bt[0] + transfer[1]*bt[1])
             + transfer[1]*(b[1] - S[1]*bt[1] + transfer[0]*bt[0])
           )
           /(1+transfer.sum());

    b[1] = (b[1] + transfer[0]*(b[0] + bt[0]) - S[1]*bt[1])
           /(1+transfer[0]);
}

// bouyancy force on faces
for(label ip = 0; ip < nParts; ip++)
{
    bf[ip] = fvc::interpolate(b[ip], "b")*mesh.Sf().component(2);
}

b.updateSum();
bf.updateSum();

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


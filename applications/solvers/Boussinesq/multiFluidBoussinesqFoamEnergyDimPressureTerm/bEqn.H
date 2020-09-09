// advective form b equation

for(label ip = 0; ip < nParts; ip++)
{
    b.ddt()[ip] = -fvc::div(volFlux[ip], bBar + b[ip], "b")
                + (bBar + b[ip])*fvc::div(volFlux[ip])
                + alpha*fvc::laplacian(b[ip]);
}

if (!adjustQ)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        b.ddt()[ip] += Q;
    }
}

if (relaxbTime.value() > 0)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        b.ddt()[ip] -= b[ip]/relaxbTime;
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
    if (transferAnomaliesLinComb)
    {
        bt[0] = bAnomalyLinearCoeff*(b[1] - b[0]);
        bt[1] = - bt[0];
    }
    
    // Bouyancy transferred is constant factor times mean buoyancy
    if (transferAnomaliesMeanFactor)
    {
        bt[0] = bAnomalyMeanCoeff*mag(b[0]);
        bt[1] = -bAnomalyMeanCoeff*mag(b[1]);
    }

    // The bouyancy of the fluid transferred depends on the laplacian of b
    if (transferAnomaliesLaplacian)
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
    // Check that bt is big enough for wTransfer technique (giving w[0]<0)
    if (wTransfer)
    {
        for(label ip = 0; ip < nParts; ip++)
        {
            // Update bMinMax
            bMinMax[ip] = -(fvc::reconstruct
            (
                volFlux.ddt()[ip]
              + volFlux[ip]*(Sf[ip]/dt + 1/wTransferTimeScale)
            ) & vector(0,0,1));
        }
        bt[0] = max
        (
            bt[0],
            (-bMinMax[0]/wTransferTimeScale + b.ddt()[0])
            /max(S[0]/dt, rateSMALL)
        );
        bt[1] = min
        (
            bt[1],
            (-bMinMax[1]/wTransferTimeScale + b.ddt()[1])
            /max(S[1]/dt, rateSMALL)
        );
    }

    // Normal transfers
    if (!transferAnomaliesb0 && !transferZeroAnomalies)
    {
        b[0] = (
                   (1+transfer[0])*(b[0] - S[0]*bt[0] + transfer[1]*bt[1])
                 + transfer[1]*(b[1] - S[1]*bt[1] + transfer[0]*bt[0])
               )
               /(1+transfer.sum());

        b[1] = (b[1] + transfer[0]*(b[0] + bt[0]) - S[1]*bt[1])
               /(1+transfer[0]);
    }
    else if (transferAnomaliesb0)
    // Fluid transferred from 1 to zero has zero buoyancy
    {
        b[0] /= 1 + transfer[1];
        //bt[0] =  remains zero
        b[1] = (b[1] + transfer[0]*b[0])/(1+transfer[0] - S[1]);
        bt[1] = -b[1];
    }
    else // transferZeroAnomalies
    // Fluid transferred has zero buoyancy
    {
        b[0] /= 1 - S[0] + transfer[1];
        b[1] /= 1 - S[1] + transfer[0];
    }
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
    dimensionedScalar totalHeat("totalHeat", dimensionSet(0, 1,-2,0,0), scalar(0));
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


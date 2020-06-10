// The transfer terms are only formulated for nParts == 2
if (nParts != 2)
{
    FatalErrorIn("massTransfer.H")
        << " mass transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}

// Set all mass transfers to zero
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=ip+1; jp < nParts; jp++)
    {
        massTransfer[ip][jp] == rateZero;
        massTransfer[jp][ip] == rateZero;
    }
}

// Mass tranfer terms as diffusion between partitions
if (Ksigma.value() > SMALL)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        for(label jp=ip+1; jp < nParts; jp++)
        {
            massTransferTmp = 0.5*Ksigma*fvc::laplacian
            (
                sigma[jp] -sigma[ip]
            );
            massTransfer[ip][jp] += max(massTransferTmp, rateZero);
            massTransfer[jp][ip] -= min(massTransferTmp, rateZero);
        }
    }
}

// Entrainment based on plume radius
if (entrainmentCoefficient > 0 && cloudRadiusMax.value() > 0)
{
    massTransfer[0][1] += max
    (
        sigma[0]*u[1].component(2)*entrainmentCoefficient/cloudRadiusMax,
        rateZero
    );
    massTransfer[1][0] += max
    (
        sigma[1]*u[0].component(2)*entrainmentCoefficient/cloudRadiusMax,
        rateZero
    );
}

// Mass transfer to avoid changes in sigma due to partition divergence
if (divTransfer>0)
{
    massTransfer[0][1] += divTransfer*sigma[0]
        *max(-fvc::div(volFlux[0]), rateZero);
    massTransfer[1][0] += divTransfer*sigma[1]
        *max(-fvc::div(volFlux[1]), rateZero);
}

// Mass transfer based on buoyancy
if (buoyancyTransfer)
{
    massTransfer[0][1] -= sigma[0]*bTransferDiffusivity*min
    (
        fvc::laplacian(b[0])/b[0],
        dimensionedScalar("", dimensionSet(0,-2,0,0,0), scalar(0))
    );

    massTransfer[1][0] += sigma[1]*bTransferDiffusivity*max
    (
        fvc::laplacian(b[1])/b[1],
        dimensionedScalar("", dimensionSet(0,-2,0,0,0), scalar(0))
    );
}

// Transfers to keep sigma uniform at selected locations
if (uniformSigma)
{
    massTransfer[0][1] = transferLocation*max(sigma.ddt()[0], rateZero)
                       + (1-transferLocation)*massTransfer[0][1];
    massTransfer[1][0] = transferLocation*max(sigma.ddt()[1], rateZero)
                       + (1-transferLocation)*massTransfer[1][0];
}

// Mass transfer based on vertical velocity and horizontal divergence
if (DwDtTransfer)
{
    DuDt[0] = fvc::reconstruct
    (
        volFlux.ddt()[0]
      + (
            mesh.Sf() & fvc::interpolate
            (
                fvc::div(volFlux[0], u[0], "u")
              - u[0]*fvc::div(volFlux[0]),
                "convection"
            )
        )
    );
    DuDt[1] = fvc::reconstruct
    (
        volFlux.ddt()[1]
      + (
            mesh.Sf() & fvc::interpolate
            (
                fvc::div(volFlux[1], u[1], "u")
              - u[1]*fvc::div(volFlux[1]),
                "convection"
            )
        )
    );

    wVar[0] = sqr(u[0].component(2) - u.mean().component(2));
    wVar[1] = sqr(u[1].component(2) - u.mean().component(2));

    massTransfer[0][1] += sigma[0]/dt * 0.5 *
    (
        1 + Foam::erf
        (
            (u[0].component(2) + dt*DuDt[0].component(2)) /
            max
            (
                sqrt(2*wVar[0]),
                dimensionedScalar("", dimensionSet(0,1,-1,0,0), SMALL)
            )
        )
    );
    massTransfer[1][0] += sigma[1]/dt * 0.5 *
    (
        1 - Foam::erf
        (
            (u[1].component(2) + dt*DuDt[1].component(2)) /
            (
                sqrt(2*wVar[1]),
                dimensionedScalar("", dimensionSet(0,1,-1,0,0), SMALL)
            )
        )
    );

    /*massTransfer[0][1] += sigma[0] * max
        (
           - DuDt0.component(2)/
             min
             (
                 u[0].component(2),
                 dimensionedScalar("", dimensionSet(0,1,-1,0,0), -SMALL)
             ),
             rateZero
        );
    massTransfer[1][0] += sigma[1] * max
        (
           - DuDt1.component(2)/
             max
             (
                 u[0].component(2),
                 dimensionedScalar("", dimensionSet(0,1,-1,0,0), SMALL)
             ),
             rateZero
        );*/
}

// Mass transfer based on movement of the fluid interface
// Requires knowledge of the fluid interface velocity
if (gradSigmaTransfer)
{
    partitionedVolVectorField uInterface = u;
    uInterface[0].component(2) = 0*uInterface[0].component(2);
    uInterface[1].component(2) = 0*uInterface[0].component(2);

    massTransfer[0][1] = max
    (
        //((u[0]-uInterface[0]) & fvc::grad(sigma[0])) - ((u[1]-uInterface[1]) & fvc::grad(sigma[1])),
        dt*((u[0] & fvc::grad(u[0])) & fvc::grad(sigma[0])),
        rateZero
    );
    massTransfer[1][0] = max
    (
        //((u[1]-uInterface[1]) & fvc::grad(sigma[1])) - ((u[0]-uInterface[0]) & fvc::grad(sigma[0])),
        dt*((u[1] & fvc::grad(u[1])) & fvc::grad(sigma[1])),
        rateZero
    );
}


// Transfers to keep w[0] <=0 and w[1] >0
if (wTransfer)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        bMinMax[ip] = -(fvc::reconstruct
        (
            volFlux.ddt()[ip]
          + volFlux[ip]*(Sf[ip]/dt + 1/wTransferTimeScale)
        ) & vector(0,0,1));
    }
    
    massTransfer[0][1] = max
    (
        massTransfer[0][1],
        sigma[0]/max(bt[0], bSMALL)*
        (
            -bMinMax[0]/wTransferTimeScale + b.ddt()[0]
        )
    );

    massTransfer[1][0] = max
    (
        massTransfer[1][0],
        sigma[1]/min(bt[1], -bSMALL)*
        (
            -bMinMax[1]/wTransferTimeScale + b.ddt()[1]
        )
    );
}


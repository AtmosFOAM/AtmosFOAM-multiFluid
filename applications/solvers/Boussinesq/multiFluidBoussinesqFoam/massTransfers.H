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

// Transfers to keep w[0] <=0 and w[1] >0
if (wTransfer)
{
    volScalarField DDtw
//        = volVectorField(-fvc::grad(P + Pi[0] - Pi.sum())).component(2)
//        + b[0];
          = volVectorField(fvc::ddt(u[0]) + (u[0] & fvc::grad(u[0])))
            .component(2);
    
    massTransfer[0][1] = max
    (
        massTransfer[0][1],
        sigma[0]/(wPrime*Foam::sqrt(2*pi))*DDtw
            *Foam::exp(-0.5*sqr(u[0].component(2)/wPrime))
    );
    
    DDtw //= volVectorField(-fvc::grad(P + Pi[1] - Pi.sum())).component(2)
         //+ b[1];
          = volVectorField(fvc::ddt(u[1]) + (u[1] & fvc::grad(u[1])))
            .component(2);

    massTransfer[1][0] = max
    (
        massTransfer[1][0],
        sigma[1]/(wPrime*Foam::sqrt(2*pi))*DDtw
            *Foam::exp(-0.5*sqr(u[1].component(2)/wPrime))
    );

    for(label ip = 0; ip < nParts; ip++)
    {
        bMinMax[ip] = (fvc::grad(P + Pi[ip] - Pi.sum()) & vector(0,0,1))
                    - u[ip].component(2)/wTransferTimeScale;
    }
    
    bt[0] = sigma[1]*min(max((bMinMax[0] - b[0]), bZero), btMax);
    bt[1] = -sigma[0]*min(max(b[1] - bMinMax[1], bZero), btMax);
}


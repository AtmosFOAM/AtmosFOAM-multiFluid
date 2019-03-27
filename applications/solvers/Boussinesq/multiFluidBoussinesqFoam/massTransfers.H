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

// Mass transfer to avoid changes in sigma due to partition divergence
if (divTransfer>0)
{
    massTransfer[0][1] += divTransfer*sigma[0]
        *max(-fvc::div(volFlux[0]), rateZero);
    massTransfer[1][0] += divTransfer*sigma[1]
        *max(-fvc::div(volFlux[1]), rateZero);
}

if (wTransfer)
{
    // Fixed transfer rate dependent on w
    massTransfer[0][1] += max(sign(u[0].component(vector::Z)), dimlessZero)
                         *sigma[0]/ transferTimeScale;
    massTransfer[1][0] += max(sign(u[1].component(vector::Z)), dimlessZero)
                         *sigma[1]/ transferTimeScale;
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


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
    massTransferD[1] = Ksigma*fvc::laplacian(sigma[0]);

    massTransferD[0] = -min(massTransferD[1], rateZero);
    massTransferD[1] -= massTransferD[0];
    
    massTransfer[0][1] += massTransferD[0];
    massTransfer[1][0] += massTransferD[1];

    massTransferD[0] /= max(sigma[0], minSigmaDiv);
    massTransferD[1] /= max(sigma[1], minSigmaDiv);
}

// Mass transfer to avoid changes in sigma due to partition divergence
if (divTransfer>0)
{
    massTransfer[0][1] += divTransfer*sigma[0]
        *max(-fvc::div(volFlux[0]), rateZero);
    massTransfer[1][0] += divTransfer*sigma[1]
        *max(-fvc::div(volFlux[1]), rateZero);
}

// Transfers to keep sigma[0] = 0 where w[0]>0 and vice verca
if (wTransfer)
{
    volScalarField w = u[0].component(2);
    volScalarField DDtw = -volVectorField
    (
        (u[0] & fvc::grad(Uf[0])) + fvc::grad(P + p[0])
    ).component(2) + b[0];
//          = volVectorField(fvc::ddt(u[0]) + (u[0] & fvc::grad(Uf[0])))
//            .component(2);

    
    forAll(DDtw, cellI)
    {
        if
        (
            (DDtw[cellI] > -w[cellI]/dt.value() && DDtw[cellI] > SMALL)
         || w[cellI] > SMALL
        )
        {
            massTransfer[0][1][cellI] += sigma[0][cellI]/dt.value();
        }
    }
    
    w = u[1].component(2);
    DDtw = -volVectorField
    (
        (u[1] & fvc::grad(Uf[1])) + fvc::grad(P + p[1])
    ).component(2) + b[1];
        //= volVectorField(fvc::ddt(u[1]) + (u[1] & fvc::grad(Uf[1])))
         //  .component(2);

    forAll(DDtw, cellI)
    {
        if 
        (
            (DDtw[cellI] < -w[cellI]/dt.value() && DDtw[cellI] < -SMALL)
          || w[cellI] < -SMALL
        )
        {
            massTransfer[1][0][cellI] += sigma[1][cellI]/dt.value();
        }
    }
}

// Transfers to keep sigma uniform at selected locations
if (uniformSigma)
{
    massTransfer[0][1] = transferLocation*max(sigma.ddt()[0], rateZero)
                       + (1-transferLocation)*massTransfer[0][1];
    massTransfer[1][0] = transferLocation*max(sigma.ddt()[1], rateZero)
                       + (1-transferLocation)*massTransfer[1][0];
}


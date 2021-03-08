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

// Mass transfer to avoid changes in sigma due to partition divergence
if (divTransfer>0)
{
    massTransfer[0][1] = divTransfer*sigma[0]
                        *max(-fvc::div(volFlux[0]), rateZero);
    massTransfer[1][0] = divTransfer*sigma[1]
                        *max(-fvc::div(volFlux[1]), rateZero);
}

if (dwdzTransfer>0)
{
    surfaceScalarField Sfkhat = mesh.Sf().component(2)/mesh.magSf();

    volScalarField dwdz = fvc::div(volFlux[0]*Sfkhat);
    massTransfer[0][1] += dwdzTransfer*sigma[0]*max(-dwdz, rateZero);

    dwdz = fvc::div(volFlux[1]*Sfkhat);
    massTransfer[1][0] = dwdzTransfer*sigma[1]*max(-dwdz, rateZero);
}

// Mass transfer based on buoyancy
if (bTransferDiffusivity.value() > SMALL)
{
    massTransfer[0][1] += -sigma[0]*bTransferDiffusivity*min
    (
        fvc::laplacian(b[0])/(b[0]+bTransferMax),
        dimensionedScalar("", dimensionSet(0,-2,0,0,0), scalar(0))
    );

    massTransfer[1][0] += sigma[1]*bTransferDiffusivity*max
    (
        fvc::laplacian(b[1])/(b[1]-bTransferMax),
        dimensionedScalar("", dimensionSet(0,-2,0,0,0), scalar(0))
    );
}

// Transfers to keep sigma[0] = 0 where w[0]>0 and vice verca
if (wTransfer > 0)
{
    volScalarField w = u[0].component(2);
    forAll(w, cellI)
    {
        if (w[cellI] > SMALL)
        {
            massTransfer[0][1][cellI] += wTransfer*sigma[0][cellI]/dt.value();
        }
    }
    
    w = u[1].component(2);
    forAll(w, cellI)
    {
        if (w[cellI] < -SMALL)
        {
            massTransfer[1][0][cellI] += wTransfer*sigma[1][cellI]/dt.value();
        }
    }
}

// Transfers to keep sigma uniform at selected locations
if (uniformSigma)
{
    massTransfer[0][1] = transferLocation
                            *max((sigma[0] - sigma[0].oldTime())/dt, rateZero)
                       + (1-transferLocation)*massTransfer[0][1];
    massTransfer[1][0] = transferLocation
                            *max((sigma[1] - sigma[1].oldTime())/dt, rateZero)
                       + (1-transferLocation)*massTransfer[1][0];
}

// Constraint to keep sigma => 0 and to ensure not more than all transferred
// in one time step (dt sigma S <= sigma)
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp = 0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            // Do not transfer out more than there exists
            massTransfer[ip][jp] = min
            (
                massTransfer[ip][jp],
                max(sigma[ip]-minSigma, dimlessZero)/dt
            );

            // Increase values that are <minSigma (shouldn't be needed)
            massTransfer[ip][jp] = max
            (
                massTransfer[ip][jp],
                max(minSigma-sigma[jp], dimlessZero)/dt
            );
        }
    }
}

// Transfer terms for other variables
transfer[0] = dt*massTransfer[0][1]/max(sigma[1], minSigmaDiv);
transfer[1] = dt*massTransfer[1][0]/max(sigma[0], minSigmaDiv);
S[0] = dt*massTransfer[0][1]/max(sigma[0], minSigmaDiv);
S[1] = dt*massTransfer[1][0]/max(sigma[1], minSigmaDiv);
transfer.updateSum();
S.updateSum();

// Transfer terms on faces
massTransferf[0] = fvc::interpolate(massTransfer[0][1], "transfer");
massTransferf[1] = fvc::interpolate(massTransfer[1][0], "transfer");

Sf[0] = dt*massTransferf[0]/max(sigmaf[0], minSigmaDiv);
Sf[1] = dt*massTransferf[1]/max(sigmaf[1], minSigmaDiv);

transferf[0] = dt*massTransferf[0]/max(sigmaf[1], minSigmaDiv);
transferf[1] = dt*massTransferf[1]/max(sigmaf[0],minSigmaDiv);
transferf.updateSum();


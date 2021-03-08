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

// Constraint to keep sigma => 0 and to ensure not more than all transferred
// in one time step (dt sigma S <= sigma)
// This is not correct for nParts > 2
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


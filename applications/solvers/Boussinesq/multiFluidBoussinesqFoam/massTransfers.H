// Set all mass transfers to zero
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=ip+1; jp < nParts; jp++)
    {
        massTransfer(ip,jp) == rateZero;
        massTransfer(jp,ip) == rateZero;
    }
}

if (transferType == divBtransfer)
{
    // This transfer term is only formulated for nParts == 3
    if (nParts != 3)
    {
        FatalErrorIn("massTransfer.H")
            << " mass transfer term divBtransfer only formulated for nParts = 3, not nParts = "
            << nParts << exit(FatalError);
    }

        // Transfer from fluid 0 to 1 or 2 if dbdz<0
        volVectorField gradb = fvc::grad(b[0]);
        volScalarField dbdz = gradb.component(2);
        volScalarField laplacianb = fvc::laplacian(b[0]);
        volScalarField transferTmp = max(sigma[0]-minSigma[0], scalar(0))
                                    *sqrt(max(-dbdz, sqr(rateZero)));
        
        forAll(dbdz, cellI)
        {
            if (dbdz[cellI] < 0)
            {
                if (laplacianb[cellI] > 0)
                {
                    massTransfer(0,1)[cellI] += transferTmp[cellI];
                }
                else
                {
                    massTransfer(0,2)[cellI] += transferTmp[cellI];
                }
            }
        }

        // Transfer out of 1 or 2 if w approaching zero
        surfaceScalarField Sfkhat = mag(mesh.Sf().component(2))/mesh.magSf();

        for(label ip = 1; ip < nParts; ip++)
        {
            volScalarField dwdz = fvc::div(volFlux[ip]*Sfkhat);
            massTransfer(ip,0) = max(sigma[ip]-minSigma[ip], scalar(0))
                                    *max(-dwdz, rateZero);
        }
}

// Mass transfer to avoid changes in sigma due to partition divergence
else if (transferType == divTransfer)
{
    if (nParts == 2)
    {
        for(label ip = 0; ip < nParts; ip++)
        {
            massTransfer(ip,1-ip) = sigma[ip]*max(-divu[ip], rateZero);
        }
    }

    else
    {
        FatalErrorIn("massTransfer.H")
            << " divTransfer not formulated for nParts = " << nParts
            << exit(FatalError);
    }
}

else if (transferType == dwdzTransfer)
{
    // This transfer term is only formulated for nParts == 2
    if (nParts != 2)
    {
        FatalErrorIn("massTransfer.H")
            << " mass transfer term dwdzTransfer only formulated for nParts = 2, not nParts = "
            << nParts << exit(FatalError);
    }

    surfaceScalarField Sfkhat = mag(mesh.Sf().component(2))/mesh.magSf();

    volScalarField dwdz = fvc::div(volFlux[0]*Sfkhat);
    massTransfer(0,1) += transferScale*max(sigma[0]-minSigma[0], scalar(0))
                        *max(-dwdz, rateZero);

    dwdz = fvc::div(volFlux[1]*Sfkhat);
    massTransfer(1,0) += transferScale*max(sigma[1]-minSigma[1], scalar(0))
                        *max(-dwdz, rateZero);
}

// Constraint to keep sigma => 0 and to ensure not more than all transferred
// in one time step (dt sigma S <= sigma)
// This is perhaps too prescriptive for nParts > 2
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp = 0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            // Do not transfer out more than exists
            massTransfer(ip,jp) = min
            (
                massTransfer(ip,jp),
                max(sigma[ip]-minSigma[ip], dimlessZero)/(dt*(nParts-1))
            );
        }
    }
}


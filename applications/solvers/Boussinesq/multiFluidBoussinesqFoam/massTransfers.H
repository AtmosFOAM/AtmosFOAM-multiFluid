// Set all mass transfers to zero
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=ip+1; jp < nParts; jp++)
    {
        massTransfer(ip,jp) == rateZero;
        massTransfer(jp,ip) == rateZero;
    }
}

if (divBtransfer)
{
    // This transfer term is only formulated for nParts == 3
    if (nParts != 3)
    {
        FatalErrorIn("massTransfer.H")
            << " mass transfer term divBtransfer only formulated for nParts = 3, not nParts = "
            << nParts << exit(FatalError);
    }

    // Transfers out of fluids 1 and 2 if w approaching zero
    surfaceScalarField k = mag(mesh.Sf().component(2))/mesh.magSf();
    surfaceScalarField horiz = k - scalar(1);
    for(label ip = 1; ip < nParts; ip++)
    {
        volScalarField transferTmp = max(-fvc::div(volFlux[ip]*k), rateZero);
        massTransfer(ip,3-ip) += max(sigma[ip]-sigmaMin[ip], scalar(0))
                                *min(transferTmp, 1/dt);
    }

/*    // Transfer out of fluid 0 to balance based on b0 anomalies
    volScalarField laplacianb = fvc::laplacian(b[0]);
    volScalarField massTransferTmp = max(-fvc::div(volFlux[0]*horiz), rateZero)
                                   + max(-fvc::div(volFlux[0]*k), rateZero);
    massTransferTmp = max(sigma[0]-sigmaMin[0], scalar(0))
                    *min(massTransferTmp, 1/dt);
    forAll(laplacianb, cellI)
    {
        if (laplacianb[cellI] > 0)
        {
            massTransfer(0,1)[cellI] = massTransfer(2,0)[cellI]
                                     + massTransferTmp[cellI];
        }
        else
        {
            massTransfer(0,2)[cellI] = massTransfer(1,0)[cellI]
                                     + massTransferTmp[cellI];
        }
    }*/
    //massTransfer(0,1) = massTransfer(2,0);
    //massTransfer(0,2) = massTransfer(1,0);
}

if (dbdzTransfer) // Doesn't really make sense. Should be removed?
{
    // This transfer term is only formulated for nParts == 2
    if (nParts == 2)
    {
        volVectorField gradb = fvc::grad(b[0]);
        volScalarField dbdz = gradb.component(2);
        massTransfer(0,1) += dbdzTransfer*max(sigma[0]-sigmaMin[0], scalar(0))
                            *sqrt(max(-dbdz, sqr(rateZero)));

        gradb = fvc::grad(b[1]);
        dbdz = gradb.component(2);
        massTransfer(1,0) += dbdzTransfer*max(sigma[1]-sigmaMin[1], scalar(0))
                        *sqrt(max(dbdz, sqr(rateZero)));
    }
    else if (nParts == 3)
    {
        // Transfer from fluid 0 to 1 or 2 if dbdz<0
        volVectorField gradb = fvc::grad(b[0]);
        volScalarField dbdz = gradb.component(2);
        volScalarField laplacianb = fvc::laplacian(b[0]);
        volScalarField transferTmp = max(sigma[0]-sigmaMin[0], scalar(0))
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
        
        // Transfer out if w approaching zero
        surfaceScalarField Sfkhat = mesh.Sf().component(2)/mesh.magSf();

        volScalarField signw = sign(u[0].component(2));
        volScalarField dwdz = max(sigma[0]-sigmaMin[0], scalar(0))
                            *fvc::div(volFlux[0]*Sfkhat);
        massTransfer(0,1) += max(-dwdz, rateZero)*max(-signw, scalar(0));
        massTransfer(0,2) += max(-dwdz, rateZero)*max(signw, scalar(0));

        for(label ip = 1; ip < nParts; ip++)
        {
            volScalarField dwdz = max(sigma[ip]-sigmaMin[ip], scalar(0))
                                *fvc::div(volFlux[ip]*Sfkhat);
            massTransfer(ip,0) += max(-dwdz, rateZero);
        }
    }
    else
    {
        FatalErrorIn("massTransfer.H")
            << " mass transfer term dwdzTransfer only formulated for nParts = 2, or 3 not nParts = "
            << nParts << exit(FatalError);
    }
}

// 

// Mass transfer to avoid changes in sigma due to partition divergence
if (divTransfer>0)
{
    if (nParts == 2)
    {
        for(label ip = 0; ip < nParts; ip++)
        {
            volScalarField transferTmp = max(-fvc::div(volFlux[ip]), rateZero);
            transferTmp = max(sigma[ip]-sigmaMin[ip], scalar(0))
                            *min(transferTmp, 1/dt);
            massTransfer(ip,1-ip) = transferTmp;
        }
    }

    else if (nParts == 3)
    {
        for(label ip = 1; ip < nParts; ip++)
        {
            volScalarField transferTmp = max(-fvc::div(volFlux[ip]), rateZero);
            massTransfer(ip,3-ip) += max(sigma[ip]-sigmaMin[ip], scalar(0))
                                *min(transferTmp, 1/dt);
        }

        volScalarField w = u[0].component(2);
        volScalarField transferTmp = max(sigma[0]-sigmaMin[0], scalar(0))
                                *max(-fvc::div(volFlux[0]), rateZero);
        forAll(w, cellI)
        {
            label jp = w[cellI] > 0 ? 1 : 2;
            massTransfer(0,jp)[cellI] += transferTmp[cellI];
        }
    }
    else
    {
        FatalErrorIn("massTransfer.H")
            << " divTransfer not formulated for nParts = " << nParts
            << exit(FatalError);
    }
}

if (dwdzTransfer>0)
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
    massTransfer(0,1) += dwdzTransfer*max(sigma[0]-sigmaMin[0], scalar(0))
                        *max(-dwdz, rateZero);

    dwdz = fvc::div(volFlux[1]*Sfkhat);
    massTransfer(1,0) += dwdzTransfer*max(sigma[1]-sigmaMin[1], scalar(0))
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
                max(sigma[ip]-sigmaMin[ip], dimlessZero)/(dt*(nParts-1))
            );
        }
    }
}


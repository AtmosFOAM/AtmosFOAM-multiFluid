// Set all mass transfers to zero
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=ip+1; jp < nParts; jp++)
    {
        massTransfer(ip,jp) == rateZero;
        massTransfer(jp,ip) == rateZero;
    }
}

// Mass transfer to avoid changes in sigma due to partition divergence
if (divTransfer>0)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        for(label jp=ip+1; jp < nParts; jp++)
        {
            volScalarField transferTmp = 
            (
                sigma[ip]*fvc::div(volFlux[ip])
              - sigma[jp]*fvc::div(volFlux[jp])
            )/nParts;
            massTransfer(ip,jp) = -min(transferTmp, rateZero);
            massTransfer(jp,ip) =  max(transferTmp, rateZero);
        }
    }
}

if (dbdzTransfer>0)
{
    // This transfer term is only formulated for nParts == 2
    if (nParts != 2)
    {
        FatalErrorIn("massTransfer.H")
            << " mass transfer term dwdzTransfer only formulated for nParts = 2, not nParts = "
            << nParts << exit(FatalError);
    }

    volVectorField gradb = fvc::grad(b[0]);
    volScalarField dbdz = gradb.component(2);
    massTransfer(0,1) += dbdzTransfer*sigma[0]*sqrt(max(-dbdz, sqr(rateZero)));

    gradb = fvc::grad(b[1]);
    dbdz = gradb.component(2);
    massTransfer(1,0) = dbdzTransfer*sigma[1]*sqrt(max(dbdz, sqr(rateZero)));
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

    surfaceScalarField Sfkhat = mesh.Sf().component(2)/mesh.magSf();

    volScalarField dwdz = fvc::div(volFlux[0]*Sfkhat);
    massTransfer(0,1) += dwdzTransfer*sigma[0]*max(-dwdz, rateZero);

    dwdz = fvc::div(volFlux[1]*Sfkhat);
    massTransfer(1,0) = dwdzTransfer*sigma[1]*max(-dwdz, rateZero);
}

// Constraint to keep sigma => 0 and to ensure not more than all transferred
// in one time step (dt sigma S <= sigma)
// This is not sufficient for nParts > 2
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
                max(sigma[ip]-minSigma, dimlessZero)/(dt*(nParts-1))
            );
        }
    }
}


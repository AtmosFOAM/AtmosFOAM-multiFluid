// Set all mass transfers to zero
if (transferType != noTransfer) for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=ip+1; jp < nParts; jp++)
    {
        massTransfer(ip,jp) == rateZero;
        massTransfer(jp,ip) == rateZero;
    }
}

if (transferType == dwdzTransfer && nParts == 2)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        volScalarField dwdz = fvc::div(volFlux[ip]*magk);
        massTransfer(ip,1-ip) = transferScale
                                *max(sigma[ip] - minSigma[ip], scalar(0))
                                *max(-dwdz, rateZero);
    }
}

for(label ip = 0; ip < nParts; ip++)
{
    divu[ip] = fvc::div(volFlux[ip]);
}
divu.updateSum();

// Mass transfer to avoid changes in sigma due to partition divergence
if (transferType == divTransfer)
{
    for(label ip = 0; ip < nParts; ip++) for(label jp = ip+1; jp < nParts; jp++)
    {
        massTransfer(ip,jp) = 0.5*
        (
            max(-sigma[ip]*divu[ip],rateZero) + max(sigma[jp]*divu[jp],rateZero)
        );
        massTransfer(jp,ip) = 0.5*
        (
            max(-sigma[jp]*divu[jp],rateZero) + max(sigma[ip]*divu[ip],rateZero)
        );
    }
}
else if (transferType == specialTransfer)
{
    const dimensionedScalar L("L", dimLength, scalar(5000));
    const dimensionedScalar T("T", dimTime, scalar(2e3));
    massTransfer(0,1) = max
    (
        sin(mesh.C().component(2)*M_PI/L)/T, rateZero
    );
    massTransfer(1,0) = max
    (
        sin(-mesh.C().component(2)*M_PI/L)/T, rateZero
    );
}

// Constraint to keep sigma => 0 and to ensure not more than all transferred
// in one time step (dt sigma S <= sigma)
// This is perhaps too prescriptive for nParts > 2
/*for(label ip = 0; ip < nParts; ip++)
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
*/

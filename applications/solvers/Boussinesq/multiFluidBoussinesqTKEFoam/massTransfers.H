/*// Set all mass transfers to zero
if (transferType != noTransfer) for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=ip+1; jp < nParts; jp++)
    {
        massTransfer(ip,jp) == rateZero;
        massTransfer(jp,ip) == rateZero;
    }
}
*/

// Start with diffusion transfers
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=ip+1; jp < nParts; jp++)
    {
        volScalarField lapI = 2*fvc::laplacian(nuTurb[ip], sigma[ip]);
        volScalarField lapJ = 2*fvc::laplacian(nuTurb[jp], sigma[jp]);
    
        massTransfer(ip,jp)
             = max(-lapI, rateZero)
              *max( lapJ, rateZero)
              /max(lapJ - lapI, rateSMALL);
        massTransfer(jp,ip)
             = max(-lapJ, rateZero)
              *max( lapI, rateZero)
              /max(lapI - lapJ, rateSMALL);
    }
}

// Add mixing transfers
if (mixingTransfer > 0) for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=ip+1; jp < nParts; jp++)
    {
        massTransfer(ip,jp) += mixingTransfer*sigma[ip]*sigma[jp]
                    *pow(TKE[ip]*TKE[jp], 0.25)/turbLength[ip];
        massTransfer(jp,ip) += mixingTransfer*sigma[ip]*sigma[jp]
                    *pow(TKE[ip]*TKE[jp], 0.25)/turbLength[jp];
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

else if (transferType == dwdzTransfer && nParts == 3)
{
    // Transfer out of fluid 0 if unstable
    volScalarField dbdz = fvc::div(bf[0]*magk);
    volScalarField transferTmp = max
    (
        (sigma[0]-minSigma[0])/(1-minSigma[0]), scalar(0)
    )*transferScale*sqrt(max(-dbdz, sqr(rateZero)));
    
    // Transfer to the fluid that has dwdz>0
    massTransfer(0,1) += positive(fvc::div(volFlux[1]*magk))*transferTmp;
    massTransfer(0,2) += positive(fvc::div(volFlux[2]*magk))*transferTmp;
    
    // Detrainment and Entrainment
    for(label ip = 1; ip < nParts; ip++)
    {
        volScalarField dwdz = fvc::div(volFlux[ip]*magk);
        // Detrainment
        massTransfer(ip,0) += detrainRate*max(-dwdz, rateZero)*
            sigma[ip];
            //max((sigma[ip]-minSigma[ip])/(1-minSigma[ip]), scalar(0));
    
        // Entrainment
        massTransfer(0,ip) += entRate*max(dwdz, rateZero)*
            min(sigma[0], sigma[ip]);
            //max((sigma[0]-minSigma[0])/(1-minSigma[0]), scalar(0));
    }
}

// Constraint to keep sigma => 0 and to ensure not more than all transferred
// in one time step (dt sigma S <= sigma)
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


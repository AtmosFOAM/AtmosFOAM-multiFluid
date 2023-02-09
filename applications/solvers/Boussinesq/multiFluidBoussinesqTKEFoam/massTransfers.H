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
                    *pow(0.5*(TKE[ip]+TKE[jp]), 0.5)/turbLength[ip];
        massTransfer(jp,ip) += mixingTransfer*sigma[ip]*sigma[jp]
                    *pow(0.5*(TKE[ip]+TKE[jp]), 0.5)/turbLength[jp];
    }
}

if (transferType == dwdzTransfer && nParts == 2)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        volScalarField dwdz = fvc::div(volFlux[ip]*magk);
        massTransfer(ip,1-ip) = entRate*sigma[ip]*max(-dwdz, rateZero);
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
    
    // Transfer to the fluid that has greater dwdz
    volScalarField oneGreater = positive
    (
        fvc::div(volFlux[1]*magk) - fvc::div(volFlux[2]*magk)
    );
    massTransfer(0,1) += oneGreater*transferTmp;
    massTransfer(0,2) += (1-oneGreater)*transferTmp;
    
    // Detrainment and Entrainment
    for(label ip = 1; ip < nParts; ip++)
    {
        volScalarField dwdz = fvc::div(volFlux[ip]*magk);
        // Detrainment
        massTransfer(ip,0) += detrainRate*max(-dwdz, rateZero)*sigma[ip];
    
        // Entrainment
        massTransfer(0,ip) += entRate*max(dwdz, rateZero)
            *min(sigma[0], sigma[ip]);
    
        // swap 1 and 2 if in the wrong direction
        volScalarField w = u[ip].component(2);
        if (ip == 1)
        {
            massTransfer(ip,3-ip) += sigma[ip]*mag(dwdz)*positive(sign(-w));
        }
        else
        {
            massTransfer(ip,3-ip) += sigma[ip]*mag(dwdz)*positive(sign(w));
        }
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


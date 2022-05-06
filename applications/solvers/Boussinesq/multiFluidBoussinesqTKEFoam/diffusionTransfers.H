// Mass, buoyancy and momentum transfers due to diffusion

// Calculate and apply diffusive transfers
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

sigma.transferMass(massTransfer, dt);
interpolate(sigmaf, sigma);


// Bounyancy of the fluid transferred

bT[0] = b[0];
bT[1] = b[1];

// Bouyancy anomalies of fluid transferred

// Bouyancy transferred is constant factor times mean buoyancy
if (bAnomalyMeanCoeff > SMALL)
{
    bT[0] += bAnomalyMeanCoeff*sigma[0]*mag(b[0]);
    bT[1] -= bAnomalyMeanCoeff*sigma[1]*mag(b[1]);
}

// Transfer the buoyancy
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++) if (ip != jp)
    {
        b[ip] += dt/max(sigma[ip], minSigmaDiv)*
        (
            massTransfer[jp][ip]*(bT[jp] - b[ip])
          - massTransfer[ip][jp]*(bT[ip] - b[ip])
        );
    }
}


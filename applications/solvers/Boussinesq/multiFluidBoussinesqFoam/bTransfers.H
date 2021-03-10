// Bouyancy anomalies of fluid transferred

for(label ip = 0; ip < nParts; ip++) for(label jp = ip+1; jp < nParts; jp++)
{
    bT(ip,jp) = b[ip];
    bT(jp,ip) = b[jp];
}

// Bouyancy transferred is constant factor times mean buoyancy
if (bAnomalyMeanCoeff > SMALL)
{
    // only formulated for nParts == 2
    if (nParts != 2)
    {
        FatalErrorIn("bTransfer.H")
            << " bAnomalyMeanCoeff only formulated for nParts = 2, not nParts = "
            << nParts << exit(FatalError);
    }

    bT(0,1) += bAnomalyMeanCoeff*sigma[0]*mag(b[0]);
    bT(1,0) -= bAnomalyMeanCoeff*sigma[1]*mag(b[1]);
}

// Transfer the buoyancy
b.transferField(massTransfer, bT, dt);


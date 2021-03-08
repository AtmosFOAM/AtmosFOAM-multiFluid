// Bounyancy of the fluid transferred
// only formulated for nParts == 2
if (nParts != 2)
{
    FatalErrorIn("massTransfer.H")
        << " mass transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}

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
b.transferField(massTransfer, bT, dt);


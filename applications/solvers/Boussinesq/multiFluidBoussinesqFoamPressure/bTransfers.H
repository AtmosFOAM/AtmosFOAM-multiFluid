// Bouyancy anomalies of fluid transferred

for(label ip = 0; ip < nParts; ip++) for(label jp = ip+1; jp < nParts; jp++)
{
    bT(ip,jp) == b[ip];
    bT(jp,ip) == b[jp];
}

// Bouyancy transferred is constant factor times mean buoyancy
if (bTransferScale > SMALL && nParts == 2)
{
    bT(0,1) += bTransferScale*sigma[0]*max(-b[0], bZero);
    bT(1,0) -= bTransferScale*sigma[1]*max(b[1], bZero);
}
/*else if (bTransferScale > SMALL && nParts == 3)
{
    bT(0,2) -= bAnomalyMeanCoeff*sigma[0]*mag(b[0]);
    bT(2,0) += bAnomalyMeanCoeff*sigma[2]*mag(b[2]);
    bT(1,2) -= bAnomalyMeanCoeff*sigma[1]*mag(b[1]);
    bT(2,1) += bAnomalyMeanCoeff*sigma[2]*mag(b[2]);
}*/

// Transfer the buoyancy
b.transferField(massTransfer, bT, dt);


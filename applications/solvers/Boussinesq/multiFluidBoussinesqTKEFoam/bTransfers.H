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
else if (transferType == dwdzTransfer && nParts == 3 && bTransferScale > SMALL)
{
    //bT(0,1) += bTransferScale*sigma[0]*max(b[0], bZero);
    //bT(0,2) += bTransferScale*sigma[0]*min(b[0], bZero);
    bT(1,0) -= bTransferScale*sigma[1]*max(b[1], bZero);
    bT(1,2) -= bTransferScale*sigma[1]*max(b[1], bZero);
    bT(2,0) += bTransferScale*sigma[2]*max(-b[2], bZero);
    bT(2,1) += bTransferScale*sigma[2]*max(-b[2], bZero);

/*    bT(0,1) += bTransferScale*sigma[0]
             *(fvc::localMax(maxInterp.interpolate(b[0])) - b[0]);
    bT(0,2) += bTransferScale*sigma[0]
             *(fvc::localMin(-maxInterp.interpolate(-b[0])) - b[0]);
    bT(1,0) += bTransferScale*sigma[1]
             *(fvc::localMin(-maxInterp.interpolate(-b[1])) - b[1]);
    bT(1,2) += bTransferScale*sigma[1]
             *(fvc::localMin(-maxInterp.interpolate(-b[1])) - b[1]);
    bT(2,0) += bTransferScale*sigma[2]
             *(fvc::localMax(maxInterp.interpolate(b[2])) - b[2]);
    bT(2,1) += bTransferScale*sigma[2]
             *(fvc::localMax(maxInterp.interpolate(b[2])) - b[2]);
*/
}

// Transfer the buoyancy
b.transferField(massTransfer, bT, dt);


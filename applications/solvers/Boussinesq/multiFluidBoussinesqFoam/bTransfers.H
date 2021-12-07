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
else if (transferType == dwdzTransfer && bTransferScale > SMALL && nParts == 3)
{
    bT(2,0) += bTransferScale*sigma[2]*max(-b[2],bZero);
    bT(1,0) -= bTransferScale*sigma[1]*max(b[1],bZero);
    bT(0,1) += bTransferScale*sigma[0]*max(b[0],bZero);
    bT(0,2) -= bTransferScale*sigma[0]*max(b[0],bZero);
}
else if (transferType == dbdzTransfer && nParts == 3 && bTransferScale > SMALL)
{
    bT(0,1) += bTransferScale*sigma[0]
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
    
/*    bT(1,2) += min(sigma[1]*b.ddt()[1], dbdtZero)/max(massTransfer(1,2), rateSMALL);
    bT(2,1) += max(sigma[2]*b.ddt()[2], dbdtZero)/max(massTransfer(2,1), rateSMALL);

    bT(1,0) += min(sigma[1]*b.ddt()[1], dbdtZero)/max(massTransfer(1,0), rateSMALL);
    bT(2,0) += max(sigma[2]*b.ddt()[2], dbdtZero)/max(massTransfer(2,0), rateSMALL);
    
    bT(0,1) += max(sigma[0]*b.ddt()[0] + massTransfer(2,0)*(bT(2,0)-b[0]), dbdtZero)
              / max(massTransfer(0,1), rateSMALL);
    bT(0,2) += min(sigma[0]*b.ddt()[0] + massTransfer(1,0)*(bT(1,0)-b[0]), dbdtZero)
              / max(massTransfer(0,2), rateSMALL);
*/
}

// Transfer the buoyancy
b.transferField(massTransfer, bT, dt);


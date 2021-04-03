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
else if (bTransferScale > SMALL && nParts == 3)
{
    bT(2,0) += bTransferScale*sigma[2]*max(-b[2],bZero);
    bT(1,0) -= bTransferScale*sigma[1]*max(b[1],bZero);
    
    volScalarField transferTmp = sigma[0]*b.ddt()[0]
                               + massTransfer(1,0)*(bT(1,0)-b[0])
                               + massTransfer(2,0)*(bT(2,0)-b[0]);
    bT(0,1) += max(transferTmp, dbdtZero)/max(massTransfer(0,1), rateSMALL);
    bT(0,2) += min(transferTmp, dbdtZero)/max(massTransfer(0,2), rateSMALL);
}

// Transfer the buoyancy
b.transferField(massTransfer, bT, dt);


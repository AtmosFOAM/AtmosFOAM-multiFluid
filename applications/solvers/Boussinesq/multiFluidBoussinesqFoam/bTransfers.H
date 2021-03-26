// Bouyancy anomalies of fluid transferred

for(label ip = 0; ip < nParts; ip++) for(label jp = ip+1; jp < nParts; jp++)
{
    bT(ip,jp) == b[ip];
    bT(jp,ip) == b[jp];
}

if (transferType == divBtransfer && nParts == 3 && bAnomalyMeanCoeff > SMALL)
{
    bT(1,0) += sqr(sigma[1])*min(b.ddt()[1],dbdtZero)/max(massTransfer(1,0), rateSMALL);
    bT(2,0) += sqr(sigma[2])*max(b.ddt()[2],dbdtZero)/max(massTransfer(2,0), rateSMALL);
    
    volScalarField transferTmp = b.ddt()[0]
                               + massTransfer(1,0)*(bT(1,0)-b[0])
                               + massTransfer(2,0)*(bT(2,0)-b[0]);
    bT(0,1) += sqr(sigma[0])*max(transferTmp, dbdtZero)
            /max(massTransfer(0,1), rateSMALL);
    bT(0,2) += sqr(sigma[0])*min(transferTmp, dbdtZero)
            /max(massTransfer(0,2), rateSMALL);
}

// Bouyancy transferred is constant factor times mean buoyancy
else if (bAnomalyMeanCoeff > SMALL)
{
    if (nParts >= 2)
    {
        bT(0,1) += bAnomalyMeanCoeff*sigma[0]*mag(b[0]);
        bT(1,0) -= bAnomalyMeanCoeff*sigma[1]*mag(b[1]);
    }
    if (nParts >= 3)
    {
        bT(0,2) -= bAnomalyMeanCoeff*sigma[0]*mag(b[0]);
        bT(2,0) += bAnomalyMeanCoeff*sigma[2]*mag(b[2]);
        bT(1,2) -= bAnomalyMeanCoeff*sigma[1]*mag(b[1]);
        bT(2,1) += bAnomalyMeanCoeff*sigma[2]*mag(b[2]);
    }
}

// Transfer the buoyancy
b.transferField(massTransfer, bT, dt);


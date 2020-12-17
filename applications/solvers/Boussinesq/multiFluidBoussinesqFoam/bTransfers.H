// Bouyancy anomalies of fluid transferred
// Bouyancy transferred is linear combination of buoyancy in each fluid
if (bAnomalyLinearCoeff > SMALL)
{
    bt[0] = bAnomalyLinearCoeff*(b[1] - b[0]);
    bt[1] = - bt[0];
}

// Bouyancy transferred is constant factor times mean buoyancy
if (bAnomalyMeanCoeff > SMALL)
{
    bt[0] = bAnomalyMeanCoeff*sigma[0]*mag(b[0]);
    bt[1] = -bAnomalyMeanCoeff*sigma[1]*mag(b[1]);
}

// Transfer buoyancy to keep w[0] and w[1] separate
if (keepWseparate)
{
    bt[0] = sigma[0]
        *max(b[0] - volVectorField(fvc::grad(P+p[0])).component(2), bZero);
    bt[1] =-sigma[1]
        *max(b[1] - volVectorField(fvc::grad(P+p[1])).component(2), bZero);
}


// Transfers
b[0] = (
           (1+transfer[0])*(b[0] - S[0]*bt[0] + transfer[1]*bt[1])
         + transfer[1]*(b[1] - S[1]*bt[1] + transfer[0]*bt[0])
       )
       /(1+transfer.sum());

b[1] = (b[1] + transfer[0]*(b[0] + bt[0]) - S[1]*bt[1])
       /(1+transfer[0]);


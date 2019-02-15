{
// Constraint to keep sigma > 0
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp = 0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            massTransfer[ip][jp] = min
            (
                massTransfer[ip][jp],
                max
                (
                    (max(sigma[ip], dimlessZero)-minSigma)/dt,
                    rateZero
                )
            );
        }
    }
}

// Transfer terms for other variables
const scalar minSigmaDiv = max(minSigma, VSMALL);
transfer[0] = dt*massTransfer[0][1]/max(sigma[1], minSigmaDiv);
transfer[1] = dt*massTransfer[1][0]/max(sigma[0], minSigmaDiv);

Stransfer[0] = dt*massTransfer[0][1]/max(sigma[0], minSigmaDiv);
Stransfer[1] = dt*massTransfer[1][0]/max(sigma[1], minSigmaDiv);

// Apply mass transfer terms (operator split) to sigma equation
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            sigma[ip] += dt*(massTransfer[jp][ip] - massTransfer[ip][jp]);
        }
    }
}

}

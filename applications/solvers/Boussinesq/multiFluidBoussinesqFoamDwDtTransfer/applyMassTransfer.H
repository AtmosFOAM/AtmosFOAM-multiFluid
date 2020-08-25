{
// Constraint to keep sigma > 0 and to ensure not more than half transferred
// in one time step (dt sigma S <= 0.5 sigma)
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp = 0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            // Increase zero values
            massTransfer[ip][jp] = max
            (
                massTransfer[ip][jp],
                sigma[ip]/dt*max(minSigma-sigma[jp], dimlessZero)
            );

            // Not more than half transferred out
            massTransfer[ip][jp] = min
            (
                massTransfer[ip][jp],
                0.5*sigma[ip]/dt
            );

            // Do not create zero values
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
const scalar minSigmaDiv = SMALL;
transfer[0] = dt*massTransfer[0][1]/max(sigma[1], minSigmaDiv);
transfer[1] = dt*massTransfer[1][0]/max(sigma[0], minSigmaDiv);
S[0] = dt*massTransfer[0][1]/max(sigma[0], minSigmaDiv);
S[1] = dt*massTransfer[1][0]/max(sigma[1], minSigmaDiv);
transfer.updateSum();
S.updateSum();

// Transfer terms on faces
massTransferf[0] = fvc::interpolate(massTransfer[0][1], "transfer");
massTransferf[1] = fvc::interpolate(massTransfer[1][0], "transfer");

massTransferf[0] = min
(
    massTransferf[0],
    max
    (
        (max(sigmaf[0], dimlessZero)-minSigma)/dt,
        rateZero
    )
);
massTransferf[1] = min
(
    massTransferf[1],
    max
    (
        (max(sigmaf[1], dimlessZero)-minSigma)/dt,
        rateZero
    )
);

Sf[0] = dt*massTransferf[0]/max(sigmaf[0], minSigmaDiv);
Sf[1] = dt*massTransferf[1]/max(sigmaf[1], minSigmaDiv);

transferf[0] = dt*massTransferf[0]/max(sigmaf[1], minSigmaDiv);
transferf[1] = dt*massTransferf[1]/max(sigmaf[0], minSigmaDiv);
transferf.updateSum();

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
    sigma[ip].correctBoundaryConditions();
}
}

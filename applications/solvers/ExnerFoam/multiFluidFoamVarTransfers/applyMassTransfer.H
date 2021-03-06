{
// Constraint to keep sigmaRho > 0
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
                    (max(sigmaRho[ip], densityZero)-minSigma*rho[ip])/dt,
                    massTransferZero
                )
            );
        }
    }
}

// Transfer terms for other variables
transfer[0] = dt*massTransfer[0][1]/max(sigmaRho[1], minSigmaDiv*rho[1]);
transfer[1] = dt*massTransfer[1][0]/max(sigmaRho[0], minSigmaDiv*rho[0]);

transferRate[0] = massTransfer[0][1]/max(sigmaRho[0], minSigmaDiv*rho[0]);
transferRate[1] = massTransfer[1][0]/max(sigmaRho[1], minSigmaDiv*rho[1]);

// Apply mass transfer terms (operator split) to sigmaRho equation
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            sigmaRho[ip] += dt*(massTransfer[jp][ip] - massTransfer[ip][jp]);
        }
    }

    // Update the face values
    sigmaRhof[ip] = fvc::interpolate(sigmaRho[ip], "sigmaRho");
}

}

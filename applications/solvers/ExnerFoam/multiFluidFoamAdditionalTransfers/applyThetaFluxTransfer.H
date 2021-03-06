// Set theta transfer terms to be either explicit or implicit based on alphaT
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            thetaTransferTerm[ip][jp] = 1/dt * transfer[ip][jp] /
            (
                1 + alphaT *
                (
                    transfer[ip][jp]
                  + transfer[jp][ip]
                )
            );
        }
    }
}

// Apply mass transfer terms (operator split) to sigmaRho equation
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            theta[ip] =
            (
                sigmaRhoPreTransfer[ip] * thetaPreTransfer[ip] + dt *
                (
                    sigmaRhoPreTransfer[jp] * thetaPreTransfer[jp] * thetaTransferTerm[jp][ip]
                  - sigmaRhoPreTransfer[ip] * thetaPreTransfer[ip] * thetaTransferTerm[ip][jp]
                )
            ) / max(sigmaRho[ip], min(SMALL, minSigma)*rho[ip]);

            theta.ddt()[ip] =
            (
                sigmaRhoPreTransfer[ip] * thetaPreTransfer.ddt()[ip] + dt *
                (
                    sigmaRhoPreTransfer[jp] * thetaPreTransfer.ddt()[jp] * thetaTransferTerm[jp][ip]
                  - sigmaRhoPreTransfer[ip] * thetaPreTransfer.ddt()[ip] * thetaTransferTerm[ip][jp]
                )
            ) / max(sigmaRho[ip], min(SMALL, minSigma)*rho[ip]);
        }
    }
}

theta.updateSum();

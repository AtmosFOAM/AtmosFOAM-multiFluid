// Set theta transfer terms to be either explicit or implicit based on alphaT
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            thetaTransferTerm[ip][jp] = 1/dt * transfer[ip][jp] *
            (
                sigmaRhoNumerator[ip] / sigmaRhoDenominator[jp]
            ) /
            (
                1 + alphaT *
                (
                    transfer[ip][jp] * 
                    (
                        sigmaRhoNumerator[ip] / sigmaRhoDenominator[jp]
                    )
                  + transfer[jp][ip] * 
                    (
                        sigmaRhoNumerator[jp] / sigmaRhoDenominator[ip]
                    )
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
            //Info << "min/max theta transfer: " << min(thetaTransferTerm[jp][ip]).value() << " , " << max(thetaTransferTerm[jp][ip]).value() << endl;
            theta[ip] += dt * thetaTransferTerm[jp][ip] *
            (
                thetaPreTransfer[jp]
              - thetaPreTransfer[ip]
            );

            theta.ddt()[ip] += dt * thetaTransferTerm[jp][ip] *
            (
                thetaPreTransfer.ddt()[jp]
              - thetaPreTransfer.ddt()[ip]
            );
        }
    }
}

theta.updateSum();

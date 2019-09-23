
const scalar minSigmaDiv = max(minSigma, VSMALL);
partitionedVolScalarField sigmaRhoPreTransfer = sigmaRho;
partitionedSurfaceScalarField sigmaRhofPreTransfer = sigmaRhof;

// Constraint to keep min(sigma[ip]) >= minSigma
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
            
            massTransfer[ip][jp] = max
            (
                massTransfer[ip][jp],
                max
                (
                    -(max(sigmaRho[jp], densityZero)-minSigma*rho[jp])/dt,
                    massTransferZero
                )
            );
        }
    }
}

for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            transfer[ip][jp] = dt*massTransfer[ip][jp]/max
            (
                sigmaRho[ip], minSigmaDiv*rho[ip]
            );
            transferf[ip][jp] = fvc::interpolate(transfer[ip][jp], "massTransfer");
        }
    }
}

// Update mass transfer terms to be either explicit or implicit depending on alphaC
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            //Info << "Transfer min: " << min(transfer[ip][jp]).value() << endl;
            //Info << "Transfer max: " << max(transfer[ip][jp]).value() << endl;
            massTransfer[ip][jp] = 1/dt * transfer[ip][jp] * sigmaRho[ip] /
            (
                1 + alphaM*(transfer[ip][jp] + transfer[jp][ip])
            );
            
            massTransferf[ip][jp] = 1/dt * transferf[ip][jp] * sigmaRhof[ip] /
            (
                1 + alphaM*(transferf[ip][jp] + transferf[jp][ip])
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
            sigmaRho[ip] += dt*(massTransfer[jp][ip] - massTransfer[ip][jp]);
            //sigmaRho[ip] += dt*(sigmaRho[jp]*transfer[jp][ip] - sigmaRho[ip]*transfer[ip][jp]);
        }
    }

    // Update the face values
    sigmaRhof[ip] = fvc::interpolate(sigmaRho[ip], "sigmaRho");
}    
sigmaRho.updateSum();
sigmaRhof.updateSum();

/*for(label ip = 0; ip < nParts; ip++)
{   
    sigma[ip] = sigmaRho[ip]/sigmaRho.sum();
    sigmaf[ip] = fvc::interpolate(sigma[ip], "sigma");
}
sigma.updateSum();
sigmaf.updateSum();*/


//Conditions for the treatment of sigmaRho in the transfer terms.
partitionedVolScalarField sigmaRhoNumerator = sigmaRho;
partitionedVolScalarField sigmaRhoDenominator = sigmaRho;
partitionedSurfaceScalarField sigmaRhofNumerator = sigmaRhof;
partitionedSurfaceScalarField sigmaRhofDenominator = sigmaRhof;
if (transferNumeratorOld)
{
    sigmaRhoNumerator = sigmaRhoPreTransfer;
}
if (transferDenominatorOld)
{
    sigmaRhoDenominator = sigmaRhoPreTransfer;
}

//Ensure denominator not zero
for(label ip = 0; ip < nParts; ip++)
{
    sigmaRhoNumerator[ip] = max(sigmaRhoNumerator[ip], min(SMALL, minSigma)*rho[ip]);
    sigmaRhoDenominator[ip] = max(sigmaRhoDenominator[ip], min(SMALL, minSigma)*rho[ip]);
    sigmaRhofNumerator[ip] = fvc::interpolate
    (
        sigmaRhoNumerator[ip], "massTransfer"
    );
    sigmaRhofDenominator[ip] = fvc::interpolate
    (
        sigmaRhoDenominator[ip], "massTransfer"
    );
}

//Transfer theta 
if (useAdvectiveThetaTransfer)
{
    #include "applyThetaTransfer.H"
}
else if (useFluxThetaTransfer)
{
    #include "applyThetaFluxTransfer.H"
}



rhoTheta = theta.sum();
for(label ip = 0; ip < nParts; ip++)
{
    // Update the density so that pressure matches in both partitions
    rho[ip] = rhoTheta/theta[ip];

    // update sigma from rho and sigmaRho
    sigma[ip] = sigmaRho[ip]/rho[ip];
}
sigma.updateSum();
for(label ip = 0; ip < nParts; ip++)
{
    sigma[ip] /= sigma.sum();
    sigmaRho[ip] = sigma[ip]*rho[ip];
}
sigma.updateSum();
rho.updateSum();
sigmaRho.updateSum();
theta.updateSum();

// Update sigma on the face and limit to between zero and one
for(label ip = 0; ip < nParts; ip++)
{
    sigmaf[ip] = fvc::interpolate(sigma[ip], "sigma");
    sigmaf[ip] = min(max(sigmaf[ip], scalar(0)), scalar(1));
    sigmaRhof[ip] = fvc::interpolate(sigmaRho[ip]);
}
sigmaf.updateSum();
sigmaRhof.updateSum();



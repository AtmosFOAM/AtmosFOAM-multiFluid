// The transfer terms are only formulated for nParts == 2
if (nParts != 2)
{
    FatalErrorIn("massTransfer.H")
        << " mass transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}

// Set all mass transfers to zero
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=ip+1; jp < nParts; jp++)
    {
        massTransfer[ip][jp] == massTransferZero;
        massTransfer[jp][ip] == massTransferZero;
    }
}

// Mass tranfer terms as diffusion between partitions
if (Ksigma.value() > SMALL)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        for(label jp=ip+1; jp < nParts; jp++)
        {
            massTransferTmp = 0.5*Ksigma*fvc::laplacian
            (
                sigmaRho[jp] -sigmaRho[ip]
            );
            massTransfer[ip][jp] += max(massTransferTmp, massTransferZero);
            massTransfer[jp][ip] -= min(massTransferTmp, massTransferZero);
        }
    }
}

// Mass transfer to avoid changes in sigmaRho due to partition divergence
if (divTransfer)
{
    massTransfer[0][1] += sigmaRho[0]*max(-fvc::div(volFlux[0]), rateZero);
    massTransfer[1][0] += sigmaRho[1]*max(-fvc::div(volFlux[1]), rateZero);
// Reversed
//    massTransfer[0][1] += sigmaRho[0]*max(-fvc::div(volFlux[1]), rateZero);
//    massTransfer[1][0] += sigmaRho[1]*max(-fvc::div(volFlux[0]), rateZero);
}

// Mass transfer based on vertical velocity and horizontal divergence
if (wTransfer)
{
    massTransfer[0][1] -= max(sign(u[0].component(vector::Z)), dimlessZero)*
    min
    (
        fvc::div(sigmaFlux[0] - sigmaRhof[0]*((Uf[0] & ghat)*ghatSf)),
        massTransferZero
    );
    
    massTransfer[1][0] -= min(sign(u[1].component(vector::Z)), dimlessZero)*
    max
    (
        fvc::div(sigmaFlux[1] - sigmaRhof[1]*((Uf[1] & ghat)*ghatSf)),
        massTransferZero
    );
}

// Mass transfer based on buoyancy
if (thetaTransfer)
{
    massTransfer[0][1] -= sigmaRho[0]*thetaTransferDiffusivity*min
    (
        fvc::laplacian(theta[0])/theta[0],
        dimensionedScalar("", dimensionSet(0,-2,0,0,0), scalar(0))
    );

    massTransfer[1][0] += sigmaRho[1]*thetaTransferDiffusivity*max
    (
        fvc::laplacian(theta[1])/theta[1],
        dimensionedScalar("", dimensionSet(0,-2,0,0,0), scalar(0))
    );
}

#include "applyMassTransfer.H"

// Update the sums
sigmaRho.updateSum();
transfer.updateSum();


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
   // Old divTransfer
    volScalarField T =-sigmaRho.ddt()[0] - (u[0] & fvc::grad(sigmaRho[0], "advectSigmaRho0"))
                     + sigmaRho.ddt()[1] + (u[1] & fvc::grad(sigmaRho[1], "advectSigmaRho1"));
    massTransfer[0][1] -= 0.5*min(T, massTransferZero);
    massTransfer[1][0] += 0.5*max(T, massTransferZero);

    // Alternative version
//    massTransfer[0][1] += sigmaRho[0]*max(-fvc::div(volFlux[0]), rateZero);
//    massTransfer[1][0] += sigmaRho[1]*max(-fvc::div(volFlux[1]), rateZero);
}

if (wTransfer)
{
    // Fixed transfer rate dependent on w
    massTransfer[0][1] += max(sign(u[0].component(vector::Z)), dimlessZero)
                         *sigmaRho[0]/ transferTimeScale;
    massTransfer[1][0] += max(sign(u[1].component(vector::Z)), dimlessZero)
                         *sigmaRho[1]/ transferTimeScale;
    
/*    // Mass transfer based on vertical velocity and horizontal divergence
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
*/
}

// Mass transfer based on buoyancy
if (thetaTransfer)
{
    massTransfer[0][1] -= sigmaRho[0]*thetaTransferDiffusivity*min
    (
        fvc::laplacian(theta[0])/max(theta[0], minThetaDiv),
        dimensionedScalar("", dimensionSet(0,-2,0,0,0), scalar(0))
    );

    massTransfer[1][0] += sigmaRho[1]*thetaTransferDiffusivity*max
    (
        fvc::laplacian(theta[1])/max(theta[1], minThetaDiv),
        dimensionedScalar("", dimensionSet(0,-2,0,0,0), scalar(0))
    );
}

Info << endl;
Info << "Pre-transfer diagnostics" << endl;
#include "printDiags.H"

#include "applyMassTransfer.H"
//#include "sigma.H"

Info << "transfer01 minmax: [" << min(transfer[0][1]).value() << ", " 
                               << max(transfer[0][1]).value() << "]" << endl;
Info << "transfer10 minmax: [" << min(transfer[1][0]).value() << ", " 
                               << max(transfer[1][0]).value() << "]" << endl;

Info << endl;
Info << "Post-transfer diagnostics" << endl;
#include "printDiags.H"

// Update the sums
transfer[0].updateSum();
transfer[1].updateSum();


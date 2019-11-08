/*----------------------------------------------------------------------------*/
//Conduct transfers according to analytic solution (1 fluid case)
/*----------------------------------------------------------------------------*/

sigmaRhoPreTransfer = sigmaRho;
massTransfer[0] = massTransferAnalytic[0];
massTransfer[1] = massTransferAnalytic[1];

#include "applyMassTransfer.H"

/*forAll(transfer[0], celli)
{
if ((transfer[0][celli] == 0) or (transfer[1][celli] == 0))
{
Info << "Transfer ratio: " << transfer[0][celli] << " " << transfer[1][celli] << " " << (1+transfer[1][celli])/(1+transfer.sum()[celli]) << endl;
}
}*/

#include "thetaTransferMean.H"
thetaTransferAnalyticMean = thetaTransferChange;
thetaTransferAnalyticMean.write();
thetaVarTransferAnalyticMean = thetaVarTransferChange;
thetaVarTransferAnalyticMean.write();
#include "thetaTransferZero.H"
thetaTransferAnalyticZero = thetaTransferChange;
thetaTransferAnalyticZero.write();
thetaVarTransferAnalyticZero = thetaVarTransferChange;
thetaVarTransferAnalyticZero.write();
#include "thetaTransferVar.H"
thetaTransferAnalyticVar = thetaTransferChange;
thetaTransferAnalyticVar.write();
thetaVarTransferAnalyticVar = thetaVarTransferChange;
thetaVarTransferAnalyticVar.write();
#include "thetaTransferMeanVar.H"
thetaTransferAnalyticMeanVar = thetaTransferChange;
thetaTransferAnalyticMeanVar.write();
thetaVarTransferAnalyticMeanVar = thetaVarTransferChange;
thetaVarTransferAnalyticMeanVar.write();

#include "wTransferMean.H"
wTransferAnalyticMean = wTransferChange;
wTransferAnalyticMean.write();
wVarTransferAnalyticMean = wVarTransferChange;
wVarTransferAnalyticMean.write();
#include "wTransferVar.H"
wTransferAnalyticVar = wTransferChange;
wTransferAnalyticVar.write();
wVarTransferAnalyticVar = wVarTransferChange;
wVarTransferAnalyticVar.write();
#include "wTransferMeanVar.H"
wTransferAnalyticMeanVar = wTransferChange;
wTransferAnalyticMeanVar.write();
wVarTransferAnalyticMeanVar = wVarTransferChange;
wVarTransferAnalyticMeanVar.write();
#include "wTransferZero.H"
wTransferAnalyticZero = wTransferChange;
wTransferAnalyticZero.write();
wVarTransferAnalyticZero = wVarTransferChange;
wVarTransferAnalyticZero.write();

wT.write();
sigmaRho = sigmaRhoPreTransfer;


/*----------------------------------------------------------------------------*/
//Conduct transfers according to the laplacian of the fluid mass/fraction
/*----------------------------------------------------------------------------*/

sigmaRhoPreTransfer = sigmaRho;
massTransfer[0] = massTransferSigmaDiffusion[0];
massTransfer[1] = massTransferSigmaDiffusion[1];

#include "applyMassTransfer.H"

#include "thetaTransferMean.H"
thetaTransferSigmaDiffMean = thetaTransferChange;
thetaTransferSigmaDiffMean.write();
thetaVarTransferSigmaDiffMean = thetaVarTransferChange;
thetaVarTransferSigmaDiffMean.write();
#include "thetaTransferVar.H"
thetaTransferSigmaDiffVar = thetaTransferChange;
thetaTransferSigmaDiffVar.write();
thetaVarTransferSigmaDiffVar = thetaVarTransferChange;
thetaVarTransferSigmaDiffVar.write();
#include "thetaTransferMeanVar.H"
thetaTransferSigmaDiffMeanVar = thetaTransferChange;
thetaTransferSigmaDiffMeanVar.write();
thetaVarTransferSigmaDiffMeanVar = thetaVarTransferChange;
thetaVarTransferSigmaDiffMeanVar.write();

#include "wTransferMean.H"
wTransferSigmaDiffMean = wTransferChange;
wTransferSigmaDiffMean.write();
wVarTransferSigmaDiffMean = wVarTransferChange;
wVarTransferSigmaDiffMean.write();
#include "wTransferVar.H"
wTransferSigmaDiffVar = wTransferChange;
wTransferSigmaDiffVar.write();
wVarTransferSigmaDiffVar = wVarTransferChange;
wVarTransferSigmaDiffVar.write();
#include "wTransferMeanVar.H"
wTransferSigmaDiffMeanVar = wTransferChange;
wTransferSigmaDiffMeanVar.write();
wVarTransferSigmaDiffMeanVar = wVarTransferChange;
wVarTransferSigmaDiffMeanVar.write();
#include "wTransferZero.H"
wTransferSigmaDiffZero = wTransferChange;
wTransferSigmaDiffZero.write();
wVarTransferSigmaDiffZero = wVarTransferChange;
wVarTransferSigmaDiffZero.write();

sigmaRho = sigmaRhoPreTransfer;


/*----------------------------------------------------------------------------*/
//Conduct transfers according to the laplacian of theta
/*----------------------------------------------------------------------------*/

sigmaRhoPreTransfer = sigmaRho;
massTransfer[0] = massTransferThetaDiffusion[0];
massTransfer[1] = massTransferThetaDiffusion[1];

#include "applyMassTransfer.H"

#include "thetaTransferMean.H"
thetaTransferThetaDiffMean = thetaTransferChange;
thetaTransferThetaDiffMean.write();
thetaVarTransferThetaDiffMean = thetaVarTransferChange;
thetaVarTransferThetaDiffMean.write();
#include "thetaTransferVar.H"
thetaTransferThetaDiffVar = thetaTransferChange;
thetaTransferThetaDiffVar.write();
thetaVarTransferThetaDiffVar = thetaVarTransferChange;
thetaVarTransferThetaDiffVar.write();
#include "thetaTransferMeanVar.H"
thetaTransferThetaDiffMeanVar = thetaTransferChange;
thetaTransferThetaDiffMeanVar.write();
thetaVarTransferThetaDiffMeanVar = thetaVarTransferChange;
thetaVarTransferThetaDiffMeanVar.write();

#include "wTransferMean.H"
wTransferThetaDiffMean = wTransferChange;
wTransferThetaDiffMean.write();
wVarTransferThetaDiffMean = wVarTransferChange;
wVarTransferThetaDiffMean.write();
#include "wTransferVar.H"
wTransferThetaDiffVar = wTransferChange;
wTransferThetaDiffVar.write();
wVarTransferThetaDiffVar = wVarTransferChange;
wVarTransferThetaDiffVar.write();
#include "wTransferMeanVar.H"
wTransferThetaDiffMeanVar = wTransferChange;
wTransferThetaDiffMeanVar.write();
wVarTransferThetaDiffMeanVar = wVarTransferChange;
wVarTransferThetaDiffMeanVar.write();
#include "wTransferZero.H"
wTransferThetaDiffZero = wTransferChange;
wTransferThetaDiffZero.write();
wVarTransferThetaDiffZero = wVarTransferChange;
wVarTransferThetaDiffZero.write();

sigmaRho = sigmaRhoPreTransfer;


/*----------------------------------------------------------------------------*/
//Conduct transfers according to the divergence of fluid velocity
/*----------------------------------------------------------------------------*/

sigmaRhoPreTransfer = sigmaRho;
massTransfer[0] = massTransferDivTransfer[0];
massTransfer[1] = massTransferDivTransfer[1];

#include "applyMassTransfer.H"

#include "thetaTransferMean.H"
thetaTransferDivMean = thetaTransferChange;
thetaTransferDivMean.write();
thetaVarTransferDivMean = thetaVarTransferChange;
thetaVarTransferDivMean.write();
#include "thetaTransferVar.H"
thetaTransferDivVar = thetaTransferChange;
thetaTransferDivVar.write();
thetaVarTransferDivVar = thetaVarTransferChange;
thetaVarTransferDivVar.write();
#include "thetaTransferMeanVar.H"
thetaTransferDivMeanVar = thetaTransferChange;
thetaTransferDivMeanVar.write();
thetaVarTransferDivMeanVar = thetaVarTransferChange;
thetaVarTransferDivMeanVar.write();
#include "wTransferMean.H"
wTransferDivMean = wTransferChange;
wTransferDivMean.write();
wVarTransferDivMean = wVarTransferChange;
wVarTransferDivMean.write();
#include "wTransferVar.H"
wTransferDivVar = wTransferChange;
wTransferDivVar.write();
wVarTransferDivVar = wVarTransferChange;
wVarTransferDivVar.write();
#include "wTransferMeanVar.H"
wTransferDivMeanVar = wTransferChange;
wTransferDivMeanVar.write();
wVarTransferDivMeanVar = wVarTransferChange;
wVarTransferDivMeanVar.write();
#include "wTransferZero.H"
wTransferDivZero = wTransferChange;
wTransferDivZero.write();
wVarTransferDivZero = wVarTransferChange;
wVarTransferDivZero.write();

sigmaRho = sigmaRhoPreTransfer;


// Set all mass transfers to zero
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=ip+1; jp < nParts; jp++)
    {
        massTransfer[ip][jp] == massTransferZero;
        massTransfer[jp][ip] == massTransferZero;
    }
}

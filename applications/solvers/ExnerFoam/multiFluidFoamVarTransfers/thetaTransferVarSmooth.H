//Prevent rounding errors coursing negative variance.
thetaVar[0] = max(thetaVar[0], 0*thetaVar[0]);
thetaVar[1] = max(thetaVar[1], 0*thetaVar[1]);

thetaPreTransfer = theta;
thetaVarPreTransfer = thetaVar;


partitionedVolScalarField thetaTransferBound = theta;
thetaTransferBound[0] -= sqrt(2*max(thetaVar[0]+thetaVar[1], thetaVarMin)) * erfinv(2*dt*transferRate[0]-1);
thetaTransferBound[1] += sqrt(2*max(thetaVar[0]+thetaVar[1], thetaVarMin)) * erfinv(2*dt*transferRate[1]-1);

partitionedVolScalarField exptheta = sigma;
exptheta[0] = Foam::exp
(
  - 0.5 * sqr(thetaTransferBound[0] - theta[0]) / 
    max(thetaVarMin, thetaVar[0]+thetaVar[1])
);
exptheta[1] = Foam::exp
(
  - 0.5 * sqr(thetaTransferBound[1] - theta[1]) / 
    max(thetaVarMin, thetaVar[0]+thetaVar[1])
);

//Transferred temperature and variance
partitionedVolScalarField thetaT = theta;
thetaT[0] = theta[0] + thetaVar[0] / sqrt(2*M_PI*max(thetaVarMin, thetaVar[0]+thetaVar[1])) *
    exptheta[0] / max(scalar(1e-10), dt*transferRate[0]);
thetaT[1] = theta[1] - thetaVar[1] / sqrt(2*M_PI*max(thetaVarMin, thetaVar[0]+thetaVar[1])) *
    exptheta[1] / max(scalar(1e-10), dt*transferRate[1]);

partitionedVolScalarField thetaVarT = thetaVar;
thetaVarT[0] = thetaVar[0] + sqr(theta[0]) + thetaVar[0]*
(
    2*theta[0]
  - (theta[0]-thetaTransferBound[0])*thetaVar[0]/max(thetaVarMin, thetaVar[0]+thetaVar[1])
) / sqrt(2*M_PI*max(thetaVarMin, thetaVar[0]+thetaVar[1])) *
    exptheta[0] / max(scalar(1e-10), dt*transferRate[0]);
thetaVarT[1] = thetaVar[0] + sqr(theta[1]) - thetaVar[1]*
(
    2*theta[1]
  - (theta[1]-thetaTransferBound[1])*thetaVar[1]/max(thetaVarMin, thetaVar[0]+thetaVar[1])
) / sqrt(2*M_PI*max(thetaVarMin, thetaVar[0]+thetaVar[1])) *
    exptheta[1] / max(scalar(1e-10), dt*transferRate[1]);

//Conduct the transfers
theta[0] *= max(sigmaRhoPreTransfer[0], sigmaRhoMin) / max(sigmaRho[0], sigmaRhoMin);
theta[0] -= dt*transferRate[0]*thetaT[0] * sigmaRhoPreTransfer[0] /
    max(sigmaRho[0], sigmaRhoMin);
theta[0] += dt*transferRate[1]*thetaT[1] * sigmaRhoPreTransfer[1] /
    max(sigmaRho[0], sigmaRhoMin);

theta[1] *= max(sigmaRhoPreTransfer[1], sigmaRhoMin) / max(sigmaRho[1], sigmaRhoMin);
theta[1] += dt*transferRate[0]*thetaT[0] * sigmaRhoPreTransfer[0] /
    max(sigmaRho[1], sigmaRhoMin);
theta[1] -= dt*transferRate[1]*thetaT[1] * sigmaRhoPreTransfer[1] /
    max(sigmaRho[1], sigmaRhoMin);

thetaVar[0] += sqr(thetaPreTransfer[0]); 
thetaVar[0] *= max(sigmaRhoPreTransfer[0], sigmaRhoMin) / max(sigmaRho[0], sigmaRhoMin);
thetaVar[0] -= dt*transferRate[0]*thetaVarT[0] * sigmaRhoPreTransfer[0] /
    max(sigmaRho[0], sigmaRhoMin);
thetaVar[0] += dt*transferRate[1]*thetaVarT[1] * sigmaRhoPreTransfer[1] /
    max(sigmaRho[0], sigmaRhoMin);
thetaVar[0] -= sqr(theta[0]);

thetaVar[1] += sqr(thetaPreTransfer[1]); 
thetaVar[1] *= max(sigmaRhoPreTransfer[1], sigmaRhoMin) / max(sigmaRho[1], sigmaRhoMin);
thetaVar[1] += dt*transferRate[0]*thetaVarT[0] * sigmaRhoPreTransfer[0] /
    max(sigmaRho[1], sigmaRhoMin);
thetaVar[1] -= dt*transferRate[1]*thetaVarT[1] * sigmaRhoPreTransfer[1] /
    max(sigmaRho[1], sigmaRhoMin);
thetaVar[1] -= sqr(theta[1]);


////////////////////////////////////////////////////////////////////////////////

//Prevent rounding errors coursing negative variance.
thetaVar[0] = max(thetaVar[0], 0*thetaVar[0]);
thetaVar[1] = max(thetaVar[1], 0*thetaVar[1]);

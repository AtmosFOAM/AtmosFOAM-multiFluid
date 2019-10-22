{

//Prevent rounding errors coursing negative variance.
thetaVar[0] = max(thetaVar[0], 0*thetaVar[0]);
thetaVar[1] = max(thetaVar[1], 0*thetaVar[1]);

thetaPreTransfer = theta;
thetaVarPreTransfer = thetaVar;

partitionedVolScalarField thetaTransferBound = theta;
thetaTransferBound[0] -= sqrt(2*max(thetaVar[0], thetaVarMin)) * erfinv(2*dt*transferRate[0]-1);
thetaTransferBound[1] += sqrt(2*max(thetaVar[1], thetaVarMin)) * erfinv(2*dt*transferRate[1]-1);

partitionedVolScalarField erfTheta = sigma;
erfTheta[0] *= 0;
erfTheta[0] += 1-2*dt*transferRate[0];
erfTheta[1] *= 0;
erfTheta[1] += 2*dt*transferRate[1]-1;

partitionedVolScalarField expTheta = sigma;
expTheta[0] = Foam::exp
(
  - 0.5 * sqr(thetaTransferBound[0] - theta[0]) / 
    max(thetaVarMin, thetaVar[0])
);
expTheta[1] = Foam::exp
(
  - 0.5 * sqr(thetaTransferBound[1] - theta[1]) / 
    max(thetaVarMin, thetaVar[1])
);

thetaT[0] = theta[0] + sqrt(0.5*thetaVar[0]/M_PI) / max(scalar(1e-10), dt*transferRate[0]) * 
    expTheta[0];
thetaT[1] = theta[1] - sqrt(0.5*thetaVar[1]/M_PI) / max(scalar(1e-10), dt*transferRate[1]) * 
    expTheta[1];
thetaT.write();
    

theta[0] *= max(sigmaRhoPreTransfer[0], minSigmaDiv*rho[0]) / 
    max(sigmaRho[0], minSigmaDiv*rho[0]);
theta[0] -= dt*transferRate[0]*thetaPreTransfer[0] * sigmaRhoPreTransfer[0] / 
    max(sigmaRho[0], minSigmaDiv*rho[0]);
theta[0] += dt*transferRate[1]*thetaPreTransfer[1] * sigmaRhoPreTransfer[1] / 
    max(sigmaRho[0], minSigmaDiv*rho[0]);
theta[0] -= sigmaRhoPreTransfer[0] * 0.5 * 
    sqrt(0.5*thetaVar[0]/M_PI) * expTheta[0] / 
    max(sigmaRho[0], minSigmaDiv*rho[0]);
theta[0] -= sigmaRhoPreTransfer[1] * 0.5 * 
    sqrt(0.5*thetaVar[1]/M_PI) * expTheta[1] / 
    max(sigmaRho[0], minSigmaDiv*rho[0]);

theta[1] *= max(sigmaRhoPreTransfer[1], minSigmaDiv*rho[1]) / 
    max(sigmaRho[1], minSigmaDiv*rho[1]);
theta[1] += dt*transferRate[0]*thetaPreTransfer[0] * sigmaRhoPreTransfer[0] /
    max(sigmaRho[1], minSigmaDiv*rho[1]);
theta[1] -= dt*transferRate[1]*thetaPreTransfer[1] * sigmaRhoPreTransfer[1] /
    max(sigmaRho[1], minSigmaDiv*rho[1]);
theta[1] += sigmaRhoPreTransfer[0] * 0.5 * 
    sqrt(0.5*thetaVar[0]/M_PI) * expTheta[0] /
    max(sigmaRho[1], minSigmaDiv*rho[1]);
theta[1] += sigmaRhoPreTransfer[1] * 0.5 * 
    sqrt(0.5*thetaVar[1]/M_PI) * expTheta[1] /
    max(sigmaRho[1], minSigmaDiv*rho[1]);

thetaVar[0] += sqr(thetaPreTransfer[0]); 
thetaVar[0] *= max(sigmaRhoPreTransfer[0], minSigmaDiv*rho[0]) / max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] -= dt*transferRate[0]*
    (thetaVarPreTransfer[0]+sqr(thetaPreTransfer[0])) * sigmaRhoPreTransfer[0] / 
    max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] += dt*transferRate[1]*
    (thetaVarPreTransfer[1]+sqr(thetaPreTransfer[1])) * sigmaRhoPreTransfer[1] / 
    max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] -= sigmaRhoPreTransfer[0] * 0.5 * 
    (thetaPreTransfer[0] + thetaTransferBound[0]) * 
    sqrt(0.5*thetaVarPreTransfer[0]/M_PI) * expTheta[0] / 
    max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] -= sigmaRhoPreTransfer[1] * 0.5 * 
    (thetaPreTransfer[1] + thetaTransferBound[1]) * 
    sqrt(0.5*thetaVarPreTransfer[1]/M_PI) * expTheta[1] / 
    max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] -= sqr(theta[0]);
thetaVar[0]  = max(thetaVar[0], 0*thetaVar[0]);

thetaVar[1] += sqr(thetaPreTransfer[1]);
thetaVar[1] *= max(sigmaRhoPreTransfer[1], minSigmaDiv*rho[1]) / max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] += dt*transferRate[0]*
    (thetaVarPreTransfer[0]+sqr(thetaPreTransfer[0])) * sigmaRhoPreTransfer[0] /
    max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] -= dt*transferRate[1]*
    (thetaVarPreTransfer[1]+sqr(thetaPreTransfer[1])) * sigmaRhoPreTransfer[1] /
    max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] += sigmaRhoPreTransfer[0] * 0.5 * 
    (thetaPreTransfer[0] + thetaTransferBound[0]) * 
    sqrt(0.5*thetaVarPreTransfer[0]/M_PI) * expTheta[0] /
    max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] += sigmaRhoPreTransfer[1] * 0.5 * 
    (thetaPreTransfer[1] + thetaTransferBound[1]) * 
    sqrt(0.5*thetaVarPreTransfer[1]/M_PI) * expTheta[1] /
    max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] -= sqr(theta[1]);
thetaVar[1]  = max(thetaVar[1], 0*thetaVar[1]);

//Prevent rounding errors coursing negative variance.
thetaVar[0] = max(thetaVar[0], 0*thetaVar[0]);
thetaVar[1] = max(thetaVar[1], 0*thetaVar[1]);

//Update diagnostics
thetaTransferChange[1] = (theta[1]-thetaPreTransfer[1])/dt;
thetaTransferChange[0] = (theta[0]-thetaPreTransfer[0])/dt;

thetaVarTransferChange[1] = (thetaVar[1]-thetaVarPreTransfer[1])/dt;
thetaVarTransferChange[0] = (thetaVar[0]-thetaVarPreTransfer[0])/dt;


//Reset prognostic variables
theta = thetaPreTransfer;
thetaVar = thetaVarPreTransfer;

}

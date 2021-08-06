{

volScalarField thetaVarTransferRate = 0*transferRate[1] + 1 / thetaVarTimescale;
thetaVarTransferRate = min(thetaVarTransferRate, 1/dt);
thetaVarTransferRate *= thetaVarLengthScale;

thetaVar[0] = max(thetaVar[0], 0*thetaVar[0]);
thetaVar[1] = max(thetaVar[1], 0*thetaVar[1]);
thetaPreTransfer = theta;
thetaVarPreTransfer = thetaVar;
thetaPreTransferTot = theta;
thetaVarPreTransferTot = thetaVar;

//theta transfer
theta[0] = ((1+transfer[0])*thetaPreTransfer[0] + transfer[1]*thetaPreTransfer[1])
           /(1+transfer.sum());
theta[1] = ((1+transfer[1])*thetaPreTransfer[1] + transfer[0]*thetaPreTransfer[0])
           /(1+transfer.sum());
//theta[1] = (theta[1] + transfer[0]*theta[0])/(1+transfer[0]);


//thetaVar transfer
thetaVar[0] = 
(
    (1 + transfer[0]) * (thetaVarPreTransfer[0] + sqr(thetaPreTransfer[0])) 
  + transfer[1] * (thetaVarPreTransfer[1] + sqr(thetaPreTransfer[1]))
) / 
(
    1 + transfer.sum()
) - sqr(theta[0]);

thetaVar[1] = 
(
    (1 + transfer[1]) * (thetaVarPreTransfer[1] + sqr(thetaPreTransfer[1]))
  + transfer[0] * (thetaVarPreTransfer[0] + sqr(thetaPreTransfer[0]))
) / 
(
    1 + transfer.sum()
) - sqr(theta[1]);


//Apply corrections due to theta variance.
thetaVar[0] = max(thetaVar[0], 0*thetaVar[0]);
thetaVar[1] = max(thetaVar[1], 0*thetaVar[1]);
thetaPreTransfer = theta;
thetaVarPreTransfer = thetaVar;

thetaTransferTerm[0][1] = thetaVarTransferRate * 0.5 *
(
    1 + Foam::erf
    (
        (theta[0] - theta[1]) / 
        max
        (
            thetaMin, 
            Foam::sqrt( 2 * (thetaVar[0] + thetaVar[1]) ) 
        )
    )
);
thetaTransferTerm[1][0] = thetaTransferTerm[0][1];

massTransferThetaVar[0][1] = sigmaRho[0]*thetaTransferTerm[0][1];
massTransferThetaVar[1][0] = sigmaRho[1]*thetaTransferTerm[1][0];


//Factor which prevents change in mass fractions.
thetaTransferTerm[0][1] *= sigmaRho[1] / 
    max(sigmaRho[0]+sigmaRho[1], minSigmaDiv*(rho[0]+rho[1]));
thetaTransferTerm[1][0] *= sigmaRho[0] / 
    max(sigmaRho[0]+sigmaRho[1], minSigmaDiv*(rho[0]+rho[1]));

thetaVarTransferTerm[0][1] = thetaVarTransferRate / Foam::sqrt
( 
    2*M_PI * max( thetaVarMin, thetaVar[0] + thetaVar[1] ) 
) * 
Foam::exp
(
    - 0.5 * sqr(theta[1] - theta[0]) / 
    max(thetaVarMin, thetaVar[0] + thetaVar[1])  
);
thetaVarTransferTerm[1][0] = thetaVarTransferTerm[0][1];

//Factor which prevents change in mass fractions.
thetaVarTransferTerm[0][1] *= sigmaRho[1] / 
    max(sigmaRho[0]+sigmaRho[1], minSigmaDiv*(rho[0]+rho[1]));
thetaVarTransferTerm[1][0] *= sigmaRho[0] / 
    max(sigmaRho[0]+sigmaRho[1], minSigmaDiv*(rho[0]+rho[1]));

//Theta transfer based on mean theta
theta[0] -= dt*thetaTransferTerm[0][1]*thetaPreTransfer[0] * 
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
theta[0] += dt*thetaTransferTerm[1][0]*thetaPreTransfer[1] *
    sigmaRho[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
theta[0] -= dt*thetaVarTransferTerm[0][1]*thetaVarPreTransfer[0] * 
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
theta[0] -= dt*thetaVarTransferTerm[1][0]*thetaVarPreTransfer[1] *
    sigmaRho[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
    
theta[1] += dt*thetaTransferTerm[0][1]*thetaPreTransfer[0] * 
    sigmaRho[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
theta[1] -= dt*thetaTransferTerm[1][0]*thetaPreTransfer[1] * 
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
theta[1] += dt*thetaVarTransferTerm[0][1]*thetaVarPreTransfer[0] * 
    sigmaRho[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
theta[1] += dt*thetaVarTransferTerm[1][0]*thetaVarPreTransfer[1] *
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    

thetaVar[0] -= sqr(theta[0]);
thetaVar[0] += sqr(thetaPreTransfer[0]) *
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] -= dt*thetaTransferTerm[0][1]*thetaVarPreTransfer[0] * 
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] += dt*thetaTransferTerm[1][0]*thetaVarPreTransfer[1] *
    sigmaRho[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] -= dt*thetaTransferTerm[0][1]*sqr(thetaPreTransfer[0]) * 
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] += dt*thetaTransferTerm[1][0]*sqr(thetaPreTransfer[1]) * 
    sigmaRho[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
    
thetaVar[0] -= 2*dt*thetaVarTransferTerm[0][1]*thetaVarPreTransfer[0] * thetaPreTransfer[0] *
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] -= 2*dt*thetaVarTransferTerm[1][0]*thetaVarPreTransfer[1] * thetaPreTransfer[1] *
    sigmaRho[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] +=   dt*thetaVarTransferTerm[0][1]*thetaVarPreTransfer[0] * 
    (thetaPreTransfer[0] - thetaPreTransfer[1]) *
    thetaVarPreTransfer[0] / max(thetaVarPreTransfer[0] + thetaVarPreTransfer[1], thetaVarMin) *
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[0] -=   dt*thetaVarTransferTerm[1][0]*thetaVarPreTransfer[1] * 
    (thetaPreTransfer[0] - thetaPreTransfer[1]) *
    thetaVarPreTransfer[1] / max(thetaVarPreTransfer[0] + thetaVarPreTransfer[1], thetaVarMin) *
    sigmaRho[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);


thetaVar[1] -= sqr(theta[1]);
thetaVar[1] += sqr(thetaPreTransfer[1]) *
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] += dt*thetaTransferTerm[0][1]*thetaVarPreTransfer[0] * 
    sigmaRho[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] -= dt*thetaTransferTerm[1][0]*thetaVarPreTransfer[1] * 
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] += dt*thetaTransferTerm[0][1]*sqr(thetaPreTransfer[0]) * 
    sigmaRho[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] -= dt*thetaTransferTerm[1][0]*sqr(thetaPreTransfer[1]) * 
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    
thetaVar[1] += 2*dt*thetaVarTransferTerm[0][1]*thetaVarPreTransfer[0] * thetaPreTransfer[0] *
    sigmaRho[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] += 2*dt*thetaVarTransferTerm[1][0]*thetaVarPreTransfer[1] * thetaPreTransfer[1] *
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] -=   dt*thetaVarTransferTerm[0][1]*thetaVarPreTransfer[0] * 
    (thetaPreTransfer[0] - thetaPreTransfer[1]) *
    thetaVarPreTransfer[0] / max(thetaVarPreTransfer[0] + thetaVarPreTransfer[1], thetaVarMin) * 
    sigmaRho[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
thetaVar[1] +=   dt*thetaVarTransferTerm[1][0]*thetaVarPreTransfer[1] * 
    (thetaPreTransfer[0] - thetaPreTransfer[1]) *
    thetaVarPreTransfer[1] / max(thetaVarPreTransfer[0] + thetaVarPreTransfer[1], thetaVarMin) * 
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);


// Update the sum
theta.updateSum();

//Prevent rounding errors coursing negative variance.
thetaVar[0] = max(thetaVar[0], 0*thetaVar[0]);
thetaVar[1] = max(thetaVar[1], 0*thetaVar[1]);

for(label ip = 0; ip < nParts; ip++)
{
    Info << "theta" << ip << " post Var transfer: [" << min(theta[ip]).value() << ", " << max(theta[ip]).value() << "]" << endl;
}
for(label ip = 0; ip < nParts; ip++)
{
    Info << "thetaVar" << ip << " post Var transfer: [" << min(thetaVar[ip]).value() << ", " << max(thetaVar[ip]).value() << "]" << endl;
}
Info << "TOTAL INTERNAL ENERGY POST TRANSFER: " << sum(theta.sum()).value() << endl;


thetaTransferChange[1] = (theta[1]-thetaPreTransferTot[1])/dt;
thetaTransferChange[0] = (theta[0]-thetaPreTransferTot[0])/dt;

thetaVarTransferChange[1] = (thetaVar[1]-thetaVarPreTransferTot[1])/dt;
thetaVarTransferChange[0] = (thetaVar[0]-thetaVarPreTransferTot[0])/dt;

theta = thetaPreTransferTot;
thetaVar = thetaVarPreTransferTot;

}

Info << "Theta transfers based on theta variance:" << endl;
// The transfer terms are only formulated for nParts == 2
if (nParts != 2)
{
    FatalErrorIn("heatTransfer.H")
        << " heat transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}


const dimensionedScalar sigmaRhoMin("",sigmaRho[0].dimensions(),scalar(1e-16));
const dimensionedScalar thetaMin("",theta[0].dimensions(),scalar(1e-16));
const dimensionedScalar thetaVarMin("",thetaVar[0].dimensions(),scalar(1e-16));
const dimensionedScalar internalEnergyMin("",theta.sum().dimensions(),scalar(1e-16));


Info << "length scale " << thetaVarLengthScale << endl;
volScalarField thetaVarTransferRate = 0*transferRate[1] + 1 / thetaVarTimescale;
if (localThetaVarTransfer)
{
    thetaVarTransferRate *= dt * transferRate[1];
}
//Limit heat transfer to prevent instability.
thetaVarTransferRate = min(thetaVarTransferRate, 1/dt);
thetaVarTransferRate *= thetaVarLengthScale;

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


thetaTransferThetaVar[0][1] = (theta[1]-thetaPreTransfer[1])/dt;
thetaTransferThetaVar[1][0] = (theta[0]-thetaPreTransfer[0])/dt;
thetaVarTransferThetaVar[0][1] = (thetaVar[1]-thetaVarPreTransfer[1])/dt;
thetaVarTransferThetaVar[1][0] = (thetaVar[0]-thetaVarPreTransfer[0])/dt;

theta[0] = thetaPreTransfer[0];
theta[1] = thetaPreTransfer[1];
thetaVar[0] = thetaVarPreTransfer[0];
thetaVar[1] = thetaVarPreTransfer[1];   

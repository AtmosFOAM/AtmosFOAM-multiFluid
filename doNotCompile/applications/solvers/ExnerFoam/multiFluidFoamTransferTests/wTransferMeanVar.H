{

for (label ip=0; ip<nParts; ip++)
{
    w[ip] = u[ip].component(vector::Z);
}
//Prevent rounding errors coursing negative variance.
wVar[0] = max(wVar[0], 0*wVar[0]);
wVar[1] = max(wVar[1], 0*wVar[1]);

wPreTransfer = w;
wVarPreTransfer = wVar;
wPreTransferTot = w;
wVarPreTransferTot = wVar;

//theta transfer
transferf[0] = fvc::interpolate(transfer[0], "massTransfer");
transferf[1] = fvc::interpolate(transfer[1], "massTransfer");
transferf.updateSum();

volFlux[0] =
(
    (1 + transferf[0])*volFlux[0] + transferf[1]*volFlux[1]
) / 
(
    1 + transferf.sum()
);

volFlux[1] = 
(
    volFlux[1] + transferf[0]*volFlux[0]
) / 
(
    1 + transferf[0]
);

for (label ip=0; ip<nParts; ip++)
{
    sigmaFlux[ip] = sigmaRhof[ip]*volFlux[ip];
    u[ip] = fvc::reconstruct(volFlux[ip]);
    Uf[ip] = fvc::interpolate(u[ip], "Uf");
    w[ip] = u[ip].component(vector::Z);
}


//thetaVar transfer
wVarPreTransfer = wVar;

wVar[0] = 
(
    (1+transfer[0]) * (wVarPreTransfer[0] + sqr(wPreTransfer[0])) 
  + transfer[1] * (wVarPreTransfer[1] + sqr(wPreTransfer[1]))
) / 
(
    1 + transfer.sum()
) - sqr(w[0]);

wVar[1] = 
(
    (1+transfer[1]) * (wVarPreTransfer[1] + sqr(wPreTransfer[1])) 
  + transfer[0] * (wVarPreTransfer[0] + sqr(wPreTransfer[0]))
) / 
(
    1 + transfer.sum()
) - sqr(w[1]);


//Apply corrections according to variances
wPreTransfer = w;
wVarPreTransfer = wVar;
wVar[0] = max(wVar[0], 0*wVar[0]);
wVar[1] = max(wVar[1], 0*wVar[1]);

volScalarField wVarTransferRate = 0*transferRate[1] + 1 / wVarTimescale;
wVarTransferRate = min(wVarTransferRate, 1/dt);

wTransferTerm[0][1] = wVarTransferRate * 0.5 *
(
    1 + Foam::erf
    (
        (w[0] - w[1]) / 
        max
        (
            wMin, 
            Foam::sqrt( 2 * (wVar[0] + wVar[1]) ) 
        )
    )
);
wTransferTerm[1][0] = wTransferTerm[0][1];

massTransferWvar[0][1] = sigmaRho[0]*wTransferTerm[0][1];
massTransferWvar[1][0] = sigmaRho[1]*wTransferTerm[1][0];

//Factor which prevents change in mass fractions.
wTransferTerm[0][1] *= sigmaRho[1] / 
    max(sigmaRho[0]+sigmaRho[1], minSigmaDiv*(rho[0]+rho[1]));
wTransferTerm[1][0] *= sigmaRho[0] / 
    max(sigmaRho[0]+sigmaRho[1], minSigmaDiv*(rho[0]+rho[1]));

wVarTransferTerm[0][1] = wVarTransferRate / Foam::sqrt
( 
    2*M_PI * max( wVarMin, wVar[0] + wVar[1] ) 
) * 
Foam::exp
(
    - 0.5 * sqr(w[1] - w[0]) / 
    max(wVarMin, wVar[0] + wVar[1])  
);
wVarTransferTerm[1][0] = wVarTransferTerm[0][1];

//Factor which prevents change in mass fractions.
wVarTransferTerm[0][1] *= sigmaRho[1] / 
    max(sigmaRho[0]+sigmaRho[1], minSigmaDiv*(rho[0]+rho[1]));
wVarTransferTerm[1][0] *= sigmaRho[0] / 
    max(sigmaRho[0]+sigmaRho[1], minSigmaDiv*(rho[0]+rho[1]));

//w transfer based on mean w
w[0] -= dt*wTransferTerm[0][1]*wPreTransfer[0] * 
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
w[0] += dt*wTransferTerm[1][0]*wPreTransfer[1] *
    sigmaRho[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
w[0] -= dt*wVarTransferTerm[0][1]*wVarPreTransfer[0] * 
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
w[0] -= dt*wVarTransferTerm[1][0]*wVarPreTransfer[1] *
    sigmaRho[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
    
w[1] += dt*wTransferTerm[0][1]*wPreTransfer[0] * 
    sigmaRho[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
w[1] -= dt*wTransferTerm[1][0]*wPreTransfer[1] * 
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
w[1] += dt*wVarTransferTerm[0][1]*wVarPreTransfer[0] * 
    sigmaRho[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
w[1] += dt*wVarTransferTerm[1][0]*wVarPreTransfer[1] *
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    

wVar[0] -= sqr(w[0]);
wVar[0] += sqr(wPreTransfer[0]) *
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] -= dt*wTransferTerm[0][1]*wVarPreTransfer[0] * 
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] += dt*wTransferTerm[1][0]*wVarPreTransfer[1] *
    sigmaRho[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] -= dt*wTransferTerm[0][1]*sqr(wPreTransfer[0]) * 
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] += dt*wTransferTerm[1][0]*sqr(wPreTransfer[1]) * 
    sigmaRho[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
    
wVar[0] -= 2*dt*wVarTransferTerm[0][1]*wVarPreTransfer[0] * wPreTransfer[0] *
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] -= 2*dt*wVarTransferTerm[1][0]*wVarPreTransfer[1] * wPreTransfer[1] *
    sigmaRho[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] +=   dt*wVarTransferTerm[0][1]*wVarPreTransfer[0] * 
    (wPreTransfer[0] - wPreTransfer[1]) *
    wVarPreTransfer[0] / max(wVarPreTransfer[0] + wVarPreTransfer[1], wVarMin) *
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] -=   dt*wVarTransferTerm[1][0]*wVarPreTransfer[1] * 
    (wPreTransfer[0] - wPreTransfer[1]) *
    wVarPreTransfer[1] / max(wVarPreTransfer[0] + wVarPreTransfer[1], wVarMin) *
    sigmaRho[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);


wVar[1] -= sqr(w[1]);
wVar[1] += sqr(wPreTransfer[1]) *
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] += dt*wTransferTerm[0][1]*wVarPreTransfer[0] * 
    sigmaRho[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] -= dt*wTransferTerm[1][0]*wVarPreTransfer[1] * 
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] += dt*wTransferTerm[0][1]*sqr(wPreTransfer[0]) * 
    sigmaRho[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] -= dt*wTransferTerm[1][0]*sqr(wPreTransfer[1]) * 
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    
wVar[1] += 2*dt*wVarTransferTerm[0][1]*wVarPreTransfer[0] * wPreTransfer[0] *
    sigmaRho[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] += 2*dt*wVarTransferTerm[1][0]*wVarPreTransfer[1] * wPreTransfer[1] *
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] -=   dt*wVarTransferTerm[0][1]*wVarPreTransfer[0] * 
    (wPreTransfer[0] - wPreTransfer[1]) *
    wVarPreTransfer[0] / max(wVarPreTransfer[0] + wVarPreTransfer[1], wVarMin) * 
    sigmaRho[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] +=   dt*wVarTransferTerm[1][0]*wVarPreTransfer[1] * 
    (wPreTransfer[0] - wPreTransfer[1]) *
    wVarPreTransfer[1] / max(wVarPreTransfer[0] + wVarPreTransfer[1], wVarMin) * 
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);

//Prevent rounding errors coursing negative variance.
wVar[0] = max(wVar[0], 0*wVar[0]);
wVar[1] = max(wVar[1], 0*wVar[1]);

//Update diagnostics
wTransferChange[1] = (w[1]-wPreTransferTot[1])/dt;
wTransferChange[0] = (w[0]-wPreTransferTot[0])/dt;

wVarTransferChange[1] = (wVar[1]-wVarPreTransferTot[1])/dt;
wVarTransferChange[0] = (wVar[0]-wVarPreTransferTot[0])/dt;


//Reset prognostic variables
w = wPreTransfer;
wVar = wVarPreTransfer;

}

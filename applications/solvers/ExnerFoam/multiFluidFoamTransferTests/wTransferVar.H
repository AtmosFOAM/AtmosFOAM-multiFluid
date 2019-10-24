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


partitionedVolScalarField wTransferBound = w;
wTransferBound[0] -= sqrt(2*max(wVar[0], wVarMin)) * erfinv(2*dt*transferRate[0]-1);
wTransferBound[1] += sqrt(2*max(wVar[1], wVarMin)) * erfinv(2*dt*transferRate[1]-1);

partitionedVolScalarField expw = sigma;
expw[0] = Foam::exp
(
  - 0.5 * sqr(wTransferBound[0] - w[0]) / 
    max(wVarMin, wVar[0])
);
expw[1] = Foam::exp
(
  - 0.5 * sqr(wTransferBound[1] - w[1]) / 
    max(wVarMin, wVar[1])
);


wT[0] = w[0] + sqrt(0.5*wVar[0]/M_PI) / max(scalar(1e-10), dt*transferRate[0]) * 
    expw[0];
wT[1] = w[1] - sqrt(0.5*wVar[1]/M_PI) / max(scalar(1e-10), dt*transferRate[1]) * 
    expw[1];

w[0] *= max(sigmaRhoPreTransfer[0], sigmaRhoMin) / max(sigmaRho[0], sigmaRhoMin);
w[0] -= dt*transferRate[0]*wPreTransfer[0] * sigmaRhoPreTransfer[0] / 
    max(sigmaRho[0], sigmaRhoMin);
w[0] += dt*transferRate[1]*wPreTransfer[1] * sigmaRhoPreTransfer[1] / 
    max(sigmaRho[0], sigmaRhoMin);
w[0] -= sigmaRhoPreTransfer[0] * 0.5 * 
    sqrt(0.5*wVar[0]/M_PI) * expw[0] / 
    max(sigmaRho[0], sigmaRhoMin);
w[0] -= sigmaRhoPreTransfer[1] * 0.5 * 
    sqrt(0.5*wVar[1]/M_PI) * expw[1] / 
    max(sigmaRho[0], sigmaRhoMin);

w[1] *= max(sigmaRhoPreTransfer[1], sigmaRhoMin) / max(sigmaRho[1], sigmaRhoMin);
w[1] += dt*transferRate[0]*wPreTransfer[0] * sigmaRhoPreTransfer[0] /
    max(sigmaRho[1], sigmaRhoMin);
w[1] -= dt*transferRate[1]*wPreTransfer[1] * sigmaRhoPreTransfer[1] /
    max(sigmaRho[1], sigmaRhoMin);
w[1] += sigmaRhoPreTransfer[0] * 0.5 * 
    sqrt(0.5*wVar[0]/M_PI) * expw[0] /
    max(sigmaRho[1], sigmaRhoMin);
w[1] += sigmaRhoPreTransfer[1] * 0.5 * 
    sqrt(0.5*wVar[1]/M_PI) * expw[1] /
    max(sigmaRho[1], sigmaRhoMin);


wVar[0] += sqr(wPreTransfer[0]); 
wVar[0] *= max(sigmaRhoPreTransfer[0], sigmaRhoMin) / max(sigmaRho[0], sigmaRhoMin);
wVar[0] -= dt*transferRate[0]*
    (wVarPreTransfer[0]+sqr(wPreTransfer[0])) * sigmaRhoPreTransfer[0] / 
    max(sigmaRho[0], sigmaRhoMin);
wVar[0] += dt*transferRate[1]*
    (wVarPreTransfer[1]+sqr(wPreTransfer[1])) * sigmaRhoPreTransfer[1] / 
    max(sigmaRho[0], sigmaRhoMin);
wVar[0] -= sigmaRhoPreTransfer[0] * 0.5 * 
    (wPreTransfer[0] + wTransferBound[0]) * 
    sqrt(0.5*wVarPreTransfer[0]/M_PI) * expw[0] / 
    max(sigmaRho[0], sigmaRhoMin);
wVar[0] -= sigmaRhoPreTransfer[1] * 0.5 * 
    (wPreTransfer[1] + wTransferBound[1]) * 
    sqrt(0.5*wVarPreTransfer[1]/M_PI) * expw[1] / 
    max(sigmaRho[0], sigmaRhoMin);
wVar[0] -= sqr(w[0]);
wVar[0]  = max(wVar[0], 0*wVar[0]);

wVar[1] += sqr(wPreTransfer[1]);
wVar[1] *= max(sigmaRhoPreTransfer[1], sigmaRhoMin) / max(sigmaRho[1], sigmaRhoMin);
wVar[1] += dt*transferRate[0]*
    (wVarPreTransfer[0]+sqr(wPreTransfer[0])) * sigmaRhoPreTransfer[0] /
    max(sigmaRho[1], sigmaRhoMin);
wVar[1] -= dt*transferRate[1]*
    (wVarPreTransfer[1]+sqr(wPreTransfer[1])) * sigmaRhoPreTransfer[1] /
    max(sigmaRho[1], sigmaRhoMin);
wVar[1] += sigmaRhoPreTransfer[0] * 0.5 * 
    (wPreTransfer[0] + wTransferBound[0]) * 
    sqrt(0.5*wVarPreTransfer[0]/M_PI) * expw[0] /
    max(sigmaRho[1], sigmaRhoMin);
wVar[1] += sigmaRhoPreTransfer[1] * 0.5 * 
    (wPreTransfer[1] + wTransferBound[1]) * 
    sqrt(0.5*wVarPreTransfer[1]/M_PI) * expw[1] /
    max(sigmaRho[1], sigmaRhoMin);
wVar[1] -= sqr(w[1]);
wVar[1]  = max(wVar[1], 0*wVar[1]);

////////////////////////////////////////////////////////////////////////////////

//Prevent rounding errors coursing negative variance.
wVar[0] = max(wVar[0], 0*wVar[0]);
wVar[1] = max(wVar[1], 0*wVar[1]);


//Update diagnostics
wTransferChange[1] = (w[1]-wPreTransfer[1])/dt;
wTransferChange[0] = (w[0]-wPreTransfer[0])/dt;

wVarTransferChange[1] = (wVar[1]-wVarPreTransfer[1])/dt;
wVarTransferChange[0] = (wVar[0]-wVarPreTransfer[0])/dt;


//Reset prognostic variables
w = wPreTransfer;
wVar = wVarPreTransfer;

}

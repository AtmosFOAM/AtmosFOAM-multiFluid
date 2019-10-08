{

for (label ip=0; ip<nParts; ip++)
{
    w[ip] = u[ip].component(vector::Z);
}
wPreTransfer = w;
wVarPreTransfer = wVar;

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


//Update diagnostics
wTransferChange[1] = (w[1]-wPreTransfer[1])/dt;
wTransferChange[0] = (w[0]-wPreTransfer[0])/dt;

wVarTransferChange[1] = (wVar[1]-wVarPreTransfer[1])/dt;
wVarTransferChange[0] = (wVar[0]-wVarPreTransfer[0])/dt;


//Reset prognostic variables
w = wPreTransfer;
wVar = wVarPreTransfer;

}

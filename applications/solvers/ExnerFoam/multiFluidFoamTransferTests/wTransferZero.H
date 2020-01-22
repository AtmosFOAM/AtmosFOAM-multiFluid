{

for (label ip=0; ip<nParts; ip++)
{
    w[ip] = u[ip].component(vector::Z);
    wT[ip] *= 0;
}
wPreTransfer = w;
wVarPreTransfer = wVar;
partitionedSurfaceScalarField volFluxPreTransfer = volFlux;

//Velocity transfer
transferf[0] = fvc::interpolate(transfer[0], "massTransfer");
transferf[1] = fvc::interpolate(transfer[1], "massTransfer");
transferf.updateSum();

partitionedSurfaceScalarField volFluxt = volFlux;
for(label ip = 0; ip < nParts; ip++)
{
    volFluxt[ip] = -volFlux[ip]*mesh.Sf().component(2)/mesh.magSf();
}

volFlux[0] =
(
    (1 + transferf[0]) *
    (
        volFluxPreTransfer[0]
      - dt*fvc::interpolate(transferRate[0], "massTransfer")*volFluxt[0]
      + transferf[1]*volFluxt[1]
    )
  + transferf[1] *
    (
        volFluxPreTransfer[1]
      - dt*fvc::interpolate(transferRate[1], "massTransfer")*volFluxt[1]
      + transferf[0]*volFluxt[0]
    )
) /
(1 + transferf.sum());

volFlux[1] =
(
    (1 + transferf[1]) *
    (
        volFluxPreTransfer[1]
      - dt*fvc::interpolate(transferRate[1], "massTransfer")*volFluxt[1]
      + transferf[0]*volFluxt[0]
    )
  + transferf[0] *
    (
        volFluxPreTransfer[0]
      - dt*fvc::interpolate(transferRate[0], "massTransfer")*volFluxt[0]
      + transferf[1]*volFluxt[1]
    )
) /
(1 + transferf.sum());

/*volFlux[1] = 
(
    volFlux[1] 
  + transferf[0]*volFlux[0]
  + transferf[0]*volFluxt[0]
  - dt*fvc::interpolate(transferRate[1], "massTransfer")*volFluxt[1]
) / 
(1 + transferf[0]);*/

for (label ip=0; ip<nParts; ip++)
{
    sigmaFlux[ip] = sigmaRhof[ip]*volFlux[ip];
    u[ip] = fvc::reconstruct(volFlux[ip]);
    Uf[ip] = fvc::interpolate(u[ip], "Uf");
    w[ip] = u[ip].component(vector::Z);
}


//wVar transfer
wVar[0] =
(
    (1 + transfer[0]) *
    (
        wVarPreTransfer[0] + sqr(wPreTransfer[0])
      - dt*transferRate[0] * (wVarPreTransfer[0] + sqr(wPreTransfer[0]))
      + transfer[1] * (wVarPreTransfer[1] + sqr(wPreTransfer[1]))
    )
  + transfer[1] *
    (
        wVarPreTransfer[1] + sqr(wPreTransfer[1])
      - dt*transferRate[1] * (wVarPreTransfer[1] + sqr(wPreTransfer[1]))
      + transfer[0] * (wVarPreTransfer[0] + sqr(wPreTransfer[0]))
    )
) /
(1 + transfer.sum());

wVar[1] =
(
    (1 + transfer[1]) *
    (
        wVarPreTransfer[1] + sqr(wPreTransfer[1])
      - dt*transferRate[1] * (wVarPreTransfer[1] + sqr(wPreTransfer[1]))
      + transfer[0] * (wVarPreTransfer[0] + sqr(wPreTransfer[0]))
    )
  + transfer[0] *
    (
        wVarPreTransfer[0] + sqr(wPreTransfer[0])
      - dt*transferRate[0] * (wVarPreTransfer[0] + sqr(wPreTransfer[0]))
      + transfer[1] * (wVarPreTransfer[1] + sqr(wPreTransfer[1]))
    )
) /
(1 + transfer.sum());


//Update diagnostics
wTransferChange[1] = (w[1]-wPreTransfer[1])/dt;
wTransferChange[0] = (w[0]-wPreTransfer[0])/dt;

wVarTransferChange[1] = (wVar[1]-wVarPreTransfer[1])/dt;
wVarTransferChange[0] = (wVar[0]-wVarPreTransfer[0])/dt;


//Reset prognostic variables
w = wPreTransfer;
wVar = wVarPreTransfer;

}

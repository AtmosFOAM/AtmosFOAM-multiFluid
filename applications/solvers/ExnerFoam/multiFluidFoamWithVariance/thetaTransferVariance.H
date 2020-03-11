//Prevent rounding errors coursing negative variance.
thetaVar[0] = max(thetaVar[0], 0*thetaVar[0]);
thetaVar[1] = max(thetaVar[1], 0*thetaVar[1]);

partitionedVolScalarField thetaVarRef = thetaVar;
if (thetaVarTransferType == "Sharp")
{
    thetaVarRef[0] *= 0;
    thetaVarRef[1] *= 0;
}
else if (thetaVarTransferType == "Smooth")
{
    thetaVarRef[0] = thetaVar[1];
    thetaVarRef[1] = thetaVar[0];
}


//The boundary for what temperatures are transferred above/below
partitionedVolScalarField thetaTransferBound = theta;
thetaTransferBound[0] -= sqrt(2*max(thetaVar[0]+thetaVarRef[0], thetaVarMin)) * 
    erfinv(2*Stransfer[0]-1);
thetaTransferBound[1] += sqrt(2*max(thetaVar[1]+thetaVarRef[1], thetaVarMin)) * 
    erfinv(2*Stransfer[1]-1);

//Calculate Gaussian-like terms needed to find transferred temperature.
partitionedVolScalarField exptheta = sigma;
exptheta[0] = Foam::exp
(
  - 0.5 * sqr(thetaTransferBound[0] - theta[0]) / 
    max(thetaVarMin, thetaVar[0]+thetaVarRef[0])
);
exptheta[1] = Foam::exp
(
  - 0.5 * sqr(thetaTransferBound[1] - theta[1]) / 
    max(thetaVarMin, thetaVar[1]+thetaVarRef[1])
);


//Transferred temperature fields
partitionedVolScalarField thetaT = theta;

thetaT[0] = theta[0] + thetaVar[0] *
    exptheta[0] / max(scalar(1e-10), Stransfer[0]) / 
    sqrt(2*M_PI*max(thetaVarMin, thetaVar[0]+thetaVarRef[0]));

thetaT[1] = theta[1] - thetaVar[1] *
    exptheta[1] / max(scalar(1e-10), Stransfer[1]) / 
    sqrt(2*M_PI*max(thetaVarMin, thetaVar[1]+thetaVarRef[1]));


//Transferred temperature variance fields
partitionedVolScalarField thetaVarT = thetaVar;

thetaVarT[0] = thetaVar[0] + sqr(theta[0]) + thetaVar[0] *
    exptheta[0] / max(scalar(1e-10), Stransfer[0]) *
    (
        2*theta[0]
      - (theta[0]-thetaTransferBound[0])*thetaVar[0] /
        max(thetaVarMin, thetaVar[0]+thetaVarRef[0])
    ) / 
    sqrt( 2*M_PI*max(thetaVarMin, thetaVar[0]+thetaVarRef[0]) );

thetaVarT[1] = thetaVar[0] + sqr(theta[1]) - thetaVar[1] *
    exptheta[1] / max(scalar(1e-10), Stransfer[1]) *
    (
        2*theta[1]
      - (theta[1]-thetaTransferBound[1])*thetaVar[1] /
        max(thetaVarMin, thetaVar[1]+thetaVarRef[1])
    ) / 
    sqrt(2*M_PI*max(thetaVarMin, thetaVar[1]+thetaVarRef[1]));



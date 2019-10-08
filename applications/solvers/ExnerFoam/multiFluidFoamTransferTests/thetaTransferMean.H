{

thetaPreTransfer = theta;
thetaVarPreTransfer = thetaVar;

//theta transfer
theta[0] = ((1+transfer[0])*theta[0] + transfer[1]*theta[1])
           /(1+transfer.sum());
theta[1] = (theta[1] + transfer[0]*theta[0])/(1+transfer[0]);


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


//Update diagnostics
thetaTransferChange[1] = (theta[1]-thetaPreTransfer[1])/dt;
thetaTransferChange[0] = (theta[0]-thetaPreTransfer[0])/dt;

thetaVarTransferChange[1] = (thetaVar[1]-thetaVarPreTransfer[1])/dt;
thetaVarTransferChange[0] = (thetaVar[0]-thetaVarPreTransfer[0])/dt;


//Reset prognostic variables
theta = thetaPreTransfer;
thetaVar = thetaVarPreTransfer;

}

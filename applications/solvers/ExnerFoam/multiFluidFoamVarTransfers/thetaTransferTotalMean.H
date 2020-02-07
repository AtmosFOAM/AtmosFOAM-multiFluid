thetaPreTransfer = theta;
thetaVarPreTransfer = thetaVar;

partitionedVolScalarField thetat = theta;

thetat[0] = theta.mean() - theta[0];
thetat[1] = theta.mean() - theta[1];

//theta transfer
theta[0] =
(
    (1+transfer[0])*(thetaPreTransfer[0] - dt*transferRate[0]*thetat[0] + transfer[1]*thetat[1])
  + transfer[1]*(thetaPreTransfer[1] - dt*transferRate[1]*thetat[1] + transfer[0]*thetat[0])
) / (1+transfer.sum());
theta[1] =
(
    (1+transfer[1])*(thetaPreTransfer[1] - dt*transferRate[1]*thetat[1] + transfer[0]*thetat[0])
  + transfer[0]*(thetaPreTransfer[0] - dt*transferRate[0]*thetat[0] + transfer[1]*thetat[1])
) / (1+transfer.sum());

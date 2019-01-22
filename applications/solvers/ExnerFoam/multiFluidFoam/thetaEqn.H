// advective form theta equation
for(label ip = 0; ip < nParts; ip++)
{
    fvScalarMatrix thetaEqn
    (
        fvm::ddt(theta[ip])
      + fvm::div(volFlux[ip], theta[ip], "theta")
      - fvm::Sp(fvc::div(volFlux[ip]), theta[ip])
      - alpha*fvm::laplacian(theta[ip])
    );
    thetaEqn.solve();
}

// Mass transfer terms (for nParts = 2 only) block implicit solution
if (nParts != 2)
{
    FatalErrorIn("thetaEqn.H")
        << " transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}

thetat[0] = max
(
    min(-thetaTransferDiffusivity*dt*fvc::laplacian(theta[0]), thetatMax),
    thetatMin
);
/*thetat[1] = max
(
    min(thetaTransferDiffusivity*dt*fvc::laplacian(theta[1]), thetatMax),
    thetatMin
);*/

theta[0] =
(
    (1+transfer[0])*(theta[0] - transfer[1]*thetat[1] - Stransfer[0]*thetat[0])
   + transfer[1]*(theta[1] + transfer[0]*thetat[0] + Stransfer[1]*thetat[1])
)
/(1+transfer.sum());

theta[1] =
(
    theta[1]
  + transfer[0]*thetat[0] + Stransfer[1]*thetat[1]
  + transfer[0]*theta[0]
)/(1+transfer[0]);

// Update the sum and the faces values
for(label ip = 0; ip < nParts; ip++)
{
    thetaf[ip] = fvc::interpolate(theta[ip]);
}
theta.updateSum();
thetaf.updateSum();


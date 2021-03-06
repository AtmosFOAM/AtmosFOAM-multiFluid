// advective form theta equation

for(label ip = 0; ip < nParts; ip++)
{
    theta.ddt()[ip] = -fvc::div(volFlux[ip], theta[ip], "theta")
                    + theta[ip]*fvc::div(volFlux[ip]);
    
    // Viscous term
    if (alpha.value() > VSMALL)
    {
        theta.ddt()[ip] += alpha * fvc::laplacian(theta[ip]);
    }
    //theta.ddt()[ip] += Ktheta * fvc::laplacian(theta[ip]);
    
    // Update theta according to the rate of change
    theta[ip] = theta[ip].oldTime() + dt*
    (
        (1-offCentre)*theta.ddt()[ip].oldTime()
      + offCentre*theta.ddt()[ip]
    );
    
    

    theta[ip].correctBoundaryConditions();
}

// Mass transfer terms (for nParts = 2 only) block implicit solution
if (nParts != 2)
{
    FatalErrorIn("thetaEqn.H")
        << " transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}

thetaPreTransfer = theta;


/*volScalarField thetaMean = theta.sum()/sigmaRho.sum();

//Total fluid mean transfer
theta[0] = (theta[0] + (transfer[1]-dt*transferRate[0])*thetaMean)
            /(1 + transfer[1] - dt*transferRate[0]);
theta[1] = (theta[1] + (transfer[0]-dt*transferRate[1])*thetaMean)
            /(1 + transfer[0] - dt*transferRate[1]);*/

//Conduct temperature transfers.
//Note that if variance based transfers are selected,
//these are conducted in thetaVarEqn.H
if ( !thetaVarTransfer )
{
    if (thetaTransferTotalMean)
    {
        #include "thetaTransferTotalMean.H"
    }
    else
    {
        theta[0] = ((1+transfer[0])*theta[0] + transfer[1]*theta[1])
                   /(1+transfer.sum());
        theta[1] = (theta[1] + transfer[0]*theta[0])/(1+transfer[0]);
    }
}



theta.updateSum();

/*for(label ip = 0; ip < nParts; ip++)
{
    Info << "transferRate" << ip << " post Var transfer: [" << min(transferRate[ip]).value() << ", " << max(transferRate[ip]).value() << "]" << endl;
}
Info << "thetaMean post Var transfer: [" << min(theta.mean()).value() << ", " << max(theta.mean()).value() << "]" << endl;
for(label ip = 0; ip < nParts; ip++)
{
    Info << "theta" << ip << " post Var transfer: [" << min(theta[ip]).value() << ", " << max(theta[ip]).value() << "]" << endl;
}
for(label ip = 0; ip < nParts; ip++)
{
    Info << "thetaTransferBound" << ip << " post Var transfer: [" << min(thetaTransferBound[ip]).value() << ", " << max(thetaTransferBound[ip]).value() << "]" << endl;
}
for(label ip = 0; ip < nParts; ip++)
{
    Info << "erfTheta" << ip << " post Var transfer: [" << min(erfTheta[ip]).value() << ", " << max(erfTheta[ip]).value() << "]" << endl;
}
for(label ip = 0; ip < nParts; ip++)
{
    Info << "expTheta" << ip << " post Var transfer: [" << min(expTheta[ip]).value() << ", " << max(expTheta[ip]).value() << "]" << endl;
}*/

for(label ip = 0; ip < nParts; ip++)
{
    thetaf[ip] = fvc::interpolate(theta[ip]);
}
thetaf.updateSum();

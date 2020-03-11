// advective form thetaVar equation

for(label ip = 0; ip < nParts; ip++)
{
    thetaVar.ddt()[ip] = -fvc::div(volFlux[ip], thetaVar[ip], "thetaVar")
                    + thetaVar[ip]*fvc::div(volFlux[ip]);
    
    // Viscous term
    if (alpha.value() > VSMALL)
    {
        thetaVar.ddt()[ip] += alpha * fvc::laplacian(thetaVar[ip]);
    }
    
    // Update thetaVar according to the rate of change
    thetaVar[ip] = thetaVar[ip].oldTime() + dt*
    (
        (1-offCentre)*thetaVar.ddt()[ip].oldTime()
      + offCentre*thetaVar.ddt()[ip]
    );

    thetaVar[ip].correctBoundaryConditions();
}

//Prevent rounding errors coursing negative variance.
thetaVar[0] = max(thetaVar[0], 0*thetaVar[0]);
thetaVar[1] = max(thetaVar[1], 0*thetaVar[1]);

// Mass transfer terms (for nParts = 2 only) block implicit solution
if (nParts != 2)
{
    FatalErrorIn("thetaVarEqn.H")
        << " transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}

thetaVart[0] = thetaVarT[0] - thetaVar[0];
thetaVart[1] = thetaVarT[1] - thetaVar[1];

thetaVar[0] =
(
    (1+transfer[0])*(thetaVar[0] - transfer[1]*thetaVart[1] - Stransfer[0]*thetaVart[0])
   + transfer[1]*(thetaVar[1] + transfer[0]*thetaVart[0] + Stransfer[1]*thetaVart[1])
)
/(1+transfer.sum());

thetaVar[1] =
(
    thetaVar[1]
  + transfer[0]*thetaVart[0] + Stransfer[1]*thetaVart[1]
  + transfer[0]*thetaVar[0]
)/(1+transfer[0]);


//Prevent rounding errors coursing negative variance.
thetaVar[0] = max(thetaVar[0], 0*thetaVar[0]);
thetaVar[1] = max(thetaVar[1], 0*thetaVar[1]);

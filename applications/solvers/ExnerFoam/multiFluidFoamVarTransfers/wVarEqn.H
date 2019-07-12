// advective form w equation

for (label ip=0; ip<nParts; ip++)
{
    w[ip] = u[ip].component(vector::Z);
}


for(label ip = 0; ip < nParts; ip++)
{
    wVar.ddt()[ip] = -fvc::div(volFlux[ip], wVar[ip], "theta")
                       + wVar[ip]*fvc::div(volFlux[ip]);
    
    // Viscous term
    if (nu.value() > VSMALL)
    {
        wVar.ddt()[ip] += nu * fvc::laplacian(wVar[ip]);
    }
    //wVar.ddt()[ip] += Kw * fvc::laplacian(wVar[ip]);
    
    if (wVarDiffusion)
    {
        wVar.ddt()[ip] += fvc::laplacian(KwVar[ip], wVar[ip]);
    }
    
    // Update w according to the rate of change
    wVar[ip] = wVar[ip].oldTime() + dt*
    (
        (1-offCentre)*wVar.ddt()[ip].oldTime()
      + offCentre*wVar.ddt()[ip]
    );
    
    //wVar[0] += 0.12 * 0.5*dt/wVarTimescale * ( (u[0]-u[1]) & (u[0]-u[1]) );
    //wVar[ip] += 0.12 * 0.5*dt/thetaVarTimescale * ( (u[0]-u[1]) & (u[0]-u[1]) );
    //wVar[ip] -= 0.12 * 4 * dt/thetaVarTimescale * wVar[ip];
    
    

    wVar[ip].correctBoundaryConditions();
}

if (wVarProduction)
{
    wVar[0] = max(wVar[0], 0*wVar[0]);
    wVar[1] = max(wVar[1], 0*wVar[1]);
    wVarPreTransfer = wVar;

    wVar[0] += dt/wVarProductionTimescale * 4*sigma[0]*sigma[1] * 
    (
        sqr
        (
            (w[1]-w[0])/wVarProductionSeparation
          - sqrt(wVarPreTransfer[1])
        )
      - wVarPreTransfer[0]
    );
    wVar[1] += dt/wVarProductionTimescale * 4*sigma[0]*sigma[1] * 
    (
        sqr
        (
            (w[1]-w[0])/wVarProductionSeparation
          - sqrt(wVarPreTransfer[0])
        )
      - wVarPreTransfer[1]
    );
}

// Mass transfer terms (for nParts = 2 only) block implicit solution
if (nParts != 2)
{
    FatalErrorIn("wEqn.H")
        << " transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}

wVar[0] = max(wVar[0], 0*wVar[0]);
wVar[1] = max(wVar[1], 0*wVar[1]);

wVarPreTransfer = wVar;

Info << "transferRate" << ": [" << min(transferRate[0]).value() << ", " << max(transferRate[0]).value() << "]" << endl;

/*wVar[0] *= sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] += sqr(wPreTransfer[0]) *
    sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] -= sqr(w[0]);

wVar[0] -= dt*transferRate[0]*wVarPreTransfer[0] * 
    sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] -= dt*transferRate[0]*sqr(wPreTransfer[0]) * 
    sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);

wVar[0] += dt*transferRate[1]*wVarPreTransfer[1] *
    sigmaRhoPreTransfer[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);

wVar[0] += dt*transferRate[1]*sqr(wPreTransfer[1]) * 
    sigmaRhoPreTransfer[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);


    
wVar[1] *= sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] -= sqr(w[1]);
wVar[1] += sqr(wPreTransfer[1]) *
    sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] += dt*transferRate[0]*wVarPreTransfer[0] * 
    sigmaRhoPreTransfer[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] -= dt*transferRate[1]*wVarPreTransfer[1] * 
    sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] += dt*transferRate[0]*sqr(wPreTransfer[0]) * 
    sigmaRhoPreTransfer[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] -= dt*transferRate[1]*sqr(wPreTransfer[1]) * 
    sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
*/

    
if (directVarianceTransfer)
{
    for(label ip = 0; ip < partNames.size(); ip++)
    {
        Info << "w" << ip << " pre transfer: [" << min(w[ip]).value() << ", " << max(w[ip]).value() << "]" << endl;
        Info << "wVar" << ip << " pre Var transfer: [" << min(wVar[ip]).value() << ", " << max(wVar[ip]).value() << "]" << endl;
    }

    wPreTransfer = w;
    wVarPreTransfer = wVar;

    ////////////////////////////////////////////////////////////////////////////////
    // Perform the transfers
    ////////////////////////////////////////////////////////////////////////////////

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

    for(label ip = 0; ip < partNames.size(); ip++)
    {
        volVectorField du = (w[ip] - wPreTransfer[ip]) * Foam::vector(0,0,1);
        volFlux[ip] += (fvc::interpolate( du, "Uf") & mesh.Sf());
        sigmaFlux[ip] = sigmaRhof[ip]*volFlux[ip];
        u[ip] = fvc::reconstruct(volFlux[ip]);
        Uf[ip] = fvc::interpolate(u[ip], "Uf");
    }

    // Update the sum
    w.updateSum();
    sigmaFlux.updateSum();
    volFlux.updateSum();
    u.updateSum();
    Uf.updateSum();
}
else
{
    wVarPreTransfer = wVar;
    wVarPreTransfer[0] += sqr(wPreTransfer[0]);
    wVarPreTransfer[1] += sqr(wPreTransfer[1]);

    wVar[0] = ((1+transfer[0])*wVarPreTransfer[0] + transfer[1]*wVarPreTransfer[1])
               /(1+transfer.sum()) - sqr(w[0]);
    wVar[1] = ((1+transfer[1])*wVarPreTransfer[1] + transfer[0]*wVarPreTransfer[0])
               /(1+transfer.sum()) - sqr(w[1]);
}


for(label ip = 0; ip < nParts; ip++)
{
    Info << "w" << ip << " post Var transfer: [" << min(w[ip]).value() << ", " << max(w[ip]).value() << "]" << endl;
}
for(label ip = 0; ip < nParts; ip++)
{
    Info << "wVar" << ip << " post Var transfer: [" << min(wVar[ip]).value() << ", " << max(wVar[ip]).value() << "]" << endl;
}
    
wVar.updateSum();



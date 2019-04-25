Info << "Velocity transfers based on vertical velocity variance:" << endl;
// The transfer terms are only formulated for nParts == 2
if (nParts != 2)
{
    FatalErrorIn("velocityTransfer.H")
        << " heat transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}


const dimensionedScalar wMin("",w[0].dimensions(),scalar(1e-16));
const dimensionedScalar wVarMin("",wVar[0].dimensions(),scalar(1e-16));

for(label ip = 0; ip < partNames.size(); ip++)
{
    w[ip] = u[ip].component(vector::Z);
    Info << "w" << ip << " pre transfer: [" << min(w[ip]).value() << ", " << max(w[ip]).value() << "]" << endl;
}


//volScalarField wVarTransferRate = 0*transferRate[1] + 1 / thetaVarTimescale;
//if (localThetaVarTransfer)
//{
//    wVarTransferRate *= dt * transferRate[1];
//}
//Limit heat transfer to prevent instability.
//wVarTransferRate = min(wVarTransferRate, 1/dt);
dimensionedScalar dragTime("", dimTime, scalar(1));
volScalarField dragRate = dragTime * dragCoeff/cloudRadiusMin * mag(w[0] - w[1]);
volScalarField wVarTransferRate = sqr(dragRate)/wVarTimescale;

wVar[0] = max(wVar[0], 0*wVar[0]);
wVar[1] = max(wVar[1], 0*wVar[1]);

wPreTransfer = w;
wVarPreTransfer = wVar;

wTransferTerm[0][1] = wVarTransferRate * 0.5 *
(
    1 + Foam::erf
    (
        (w[0] - w[1]) / 
        max
        (
            wMin, 
            Foam::sqrt( 2 * (wVar[0] + wVar[1]) ) 
        )
    )
);
wTransferTerm[1][0] = wTransferTerm[0][1];

massTransferWvar[0][1] = sigmaRho[0]*wTransferTerm[0][1];
massTransferWvar[1][0] = sigmaRho[1]*wTransferTerm[1][0];

//Factor which prevents change in mass fractions.
wTransferTerm[0][1] *= sigmaRho[1] / 
    max(sigmaRho[0]+sigmaRho[1], minSigmaDiv*(rho[0]+rho[1]));
wTransferTerm[1][0] *= sigmaRho[0] / 
    max(sigmaRho[0]+sigmaRho[1], minSigmaDiv*(rho[0]+rho[1]));

wVarTransferTerm[0][1] = wVarTransferRate / Foam::sqrt
( 
    2*M_PI * max( wVarMin, wVar[0] + wVar[1] ) 
) * 
Foam::exp
(
    - 0.5 * sqr(w[1] - w[0]) / 
    max(wVarMin, wVar[0] + wVar[1])  
);
wVarTransferTerm[1][0] = wVarTransferTerm[0][1];

//Factor which prevents change in mass fractions.
wVarTransferTerm[0][1] *= sigmaRho[1] / 
    max(sigmaRho[0]+sigmaRho[1], minSigmaDiv*(rho[0]+rho[1]));
wVarTransferTerm[1][0] *= sigmaRho[0] / 
    max(sigmaRho[0]+sigmaRho[1], minSigmaDiv*(rho[0]+rho[1]));

//w transfer based on mean w
w[0] -= dt*wTransferTerm[0][1]*wPreTransfer[0] * 
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
w[0] += dt*wTransferTerm[1][0]*wPreTransfer[1] *
    sigmaRho[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
w[0] -= dt*wVarTransferTerm[0][1]*wVarPreTransfer[0] * 
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
w[0] -= dt*wVarTransferTerm[1][0]*wVarPreTransfer[1] *
    sigmaRho[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
    
w[1] += dt*wTransferTerm[0][1]*wPreTransfer[0] * 
    sigmaRho[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
w[1] -= dt*wTransferTerm[1][0]*wPreTransfer[1] * 
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
w[1] += dt*wVarTransferTerm[0][1]*wVarPreTransfer[0] * 
    sigmaRho[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
w[1] += dt*wVarTransferTerm[1][0]*wVarPreTransfer[1] *
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    

wVar[0] -= sqr(w[0]);
wVar[0] += sqr(wPreTransfer[0]) *
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] -= dt*wTransferTerm[0][1]*wVarPreTransfer[0] * 
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] += dt*wTransferTerm[1][0]*wVarPreTransfer[1] *
    sigmaRho[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] -= dt*wTransferTerm[0][1]*sqr(wPreTransfer[0]) * 
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] += dt*wTransferTerm[1][0]*sqr(wPreTransfer[1]) * 
    sigmaRho[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
    
wVar[0] -= 2*dt*wVarTransferTerm[0][1]*wVarPreTransfer[0] * wPreTransfer[0] *
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] -= 2*dt*wVarTransferTerm[1][0]*wVarPreTransfer[1] * wPreTransfer[1] *
    sigmaRho[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] +=   dt*wVarTransferTerm[0][1]*wVarPreTransfer[0] * 
    (wPreTransfer[0] - wPreTransfer[1]) *
    wVarPreTransfer[0] / max(wVarPreTransfer[0] + wVarPreTransfer[1], wVarMin) *
    sigmaRho[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
wVar[0] -=   dt*wVarTransferTerm[1][0]*wVarPreTransfer[1] * 
    (wPreTransfer[0] - wPreTransfer[1]) *
    wVarPreTransfer[1] / max(wVarPreTransfer[0] + wVarPreTransfer[1], wVarMin) *
    sigmaRho[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);


wVar[1] -= sqr(w[1]);
wVar[1] += sqr(wPreTransfer[1]) *
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] += dt*wTransferTerm[0][1]*wVarPreTransfer[0] * 
    sigmaRho[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] -= dt*wTransferTerm[1][0]*wVarPreTransfer[1] * 
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] += dt*wTransferTerm[0][1]*sqr(wPreTransfer[0]) * 
    sigmaRho[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] -= dt*wTransferTerm[1][0]*sqr(wPreTransfer[1]) * 
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    
wVar[1] += 2*dt*wVarTransferTerm[0][1]*wVarPreTransfer[0] * wPreTransfer[0] *
    sigmaRho[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] += 2*dt*wVarTransferTerm[1][0]*wVarPreTransfer[1] * wPreTransfer[1] *
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] -=   dt*wVarTransferTerm[0][1]*wVarPreTransfer[0] * 
    (wPreTransfer[0] - wPreTransfer[1]) *
    wVarPreTransfer[0] / max(wVarPreTransfer[0] + wVarPreTransfer[1], wVarMin) * 
    sigmaRho[0] / max(sigmaRho[1], minSigmaDiv*rho[1]);
wVar[1] +=   dt*wVarTransferTerm[1][0]*wVarPreTransfer[1] * 
    (wPreTransfer[0] - wPreTransfer[1]) *
    wVarPreTransfer[1] / max(wVarPreTransfer[0] + wVarPreTransfer[1], wVarMin) * 
    sigmaRho[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);


// Update the sum
w.updateSum();

//Prevent rounding errors coursing negative variance.
wVar[0] = max(wVar[0], 0*wVar[0]);
wVar[1] = max(wVar[1], 0*wVar[1]);

/*for(label ip = 0; ip < partNames.size(); ip++)
{
    volVectorField du = (w[ip] - wPreTransfer[ip]) * Foam::vector(0,0,1);
    volFlux[ip] += (fvc::interpolate( du, "Uf") & mesh.Sf());
    sigmaFlux[ip] = sigmaRhof[ip]*volFlux[ip];
    u[ip] = fvc::reconstruct(volFlux[ip]);
    Uf[ip] = fvc::interpolate(u[ip], "Uf");
}*/

sigmaFlux.updateSum();
volFlux.updateSum();
u.updateSum();
Uf.updateSum();

for(label ip = 0; ip < nParts; ip++)
{
    Info << "w" << ip << " post Var transfer: [" << min(w[ip]).value() << ", " << max(w[ip]).value() << "]" << endl;
}
for(label ip = 0; ip < nParts; ip++)
{
    Info << "wVar" << ip << " post Var transfer: [" << min(wVar[ip]).value() << ", " << max(wVar[ip]).value() << "]" << endl;
}
Info << "TOTAL INTERNAL ENERGY POST TRANSFER: " << sum(w.sum()).value() << endl;


wTransferWvar[0][1] = (w[1]-wPreTransfer[1])/dt;
wTransferWvar[1][0] = (w[0]-wPreTransfer[0])/dt;
wVarTransferWvar[0][1] = (wVar[1]-wVarPreTransfer[1])/dt;
wVarTransferWvar[1][0] = (wVar[0]-wVarPreTransfer[0])/dt;

w[0] = wPreTransfer[0];
w[1] = wPreTransfer[1];
wVar[0] = wVarPreTransfer[0];
wVar[1] = wVarPreTransfer[1];   



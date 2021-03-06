{
// Energy calculations

scalarList energy(nParts, scalar(0));
scalarList KE(nParts, scalar(0));
scalarList PE(nParts, scalar(0));
scalarList IE(nParts, scalar(0));
totalFluxPreTransfer *= 0;
totalMomentumPreTransfer *= 0;
energyFieldPreTransfer *= 0;

// Caclulate the initial energy in each partition
for(label ip = 0; ip < partNames.size(); ip++)
{
    volScalarField energyTmp = 0.5*sigmaRho[ip]*magSqr(u[ip]);
    energyFieldPreTransfer += energyTmp;
    //KE[ip] = (fvc::domainIntegrate(energyTmp)/Vtot).value();
    
    //surfaceScalarField energyTmpf = 0.5*sigmaRhof[ip]*magSqr(Uf[ip]);
    surfaceScalarField energyTmpf = 0.5*sigmaRhof[ip]*sqr(volFlux[ip]/mesh.magSf());
    //KE[ip] = 0.5*sum(energyTmpf).value();
    KE[ip] = sum(energyTmpf).value()/sigmaRhof[ip].size();
    
    energyTmp = - sigmaRho[ip]*(g & mesh.C());
    energyFieldPreTransfer += energyTmp;
    PE[ip] = (fvc::domainIntegrate(energyTmp)/Vtot).value();
    
    energyTmp = sigmaRho[ip]*theta[ip]*Exner*Cv;;
    energyFieldPreTransfer += energyTmp;
    IE[ip] = (fvc::domainIntegrate(energyTmp)/Vtot).value();
    
    energy[ip] = KE[ip] + PE[ip] + IE[ip];
    
    totalFluxPreTransfer += linearInterpolate(sigmaRho[ip])*volFlux[ip];
    totalMomentumPreTransfer += sigmaRho[ip]*u[ip];
}

KEsumPreTransfer = sum(KE);
PEsumPreTransfer = sum(PE);
IEsumPreTransfer = sum(IE);
energyTotPreTransfer = sum(energy);
}

// Read in prognostic variables
Info<< "Reading field Exner\n" << endl;
volScalarField Exner
(
    IOobject
    (
        "Exner",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField ExnerPrime
(
    IOobject("ExnerPrime", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("zero", dimless, scalar(0)),
    "zeroGradient"
);

// The partition names
const label nParts = 2;
wordList partNames(nParts);
partNames[0] = "stable";
partNames[1] = "buoyant";

Info << "Reading in sigma in each partition\n" << endl;
partitionedVolScalarField sigma
(
    "sigma", partNames, mesh, runTime.timeName()
);

Info << "Reading in rho in each partition\n" << endl;
partitionedVolScalarField rho
(
    "rho", partNames, mesh, runTime.timeName(), sigma
);
partitionedVolScalarField sigmaRho = rho.timesSigma();


Info << "Reading in theta in each partition\n" << endl;
partitionedVolScalarField theta
(
    "theta", partNames, mesh, runTime.timeName(), sigmaRho
);



Info << "Reading in thetaVar in each partition\n" << endl;
partitionedVolScalarField thetaVar
(
    "thetaVar", partNames, mesh, runTime.timeName(), sigmaRho
);

Info << "Calculating density and density in each partition and density weighted theta sum\n" << endl;
//for(label ip = 0; ip < partNames.size(); ip++)
//{
//    rho[ip] = pRef/(R*theta[ip])*pow(Exner, (1-kappa)/kappa);
//}
partitionedVolScalarField sigmaRhoInitial
(
    "sigmaRhoInitial", partNames, mesh, runTime.timeName()
);
partitionedVolScalarField sigmaRhoExpected
(
    "sigmaRhoExpected", partNames, mesh, runTime.timeName()
);

for(label ip = 0; ip < partNames.size(); ip++)
{
    sigmaRho[ip] = sigmaRhoInitial[ip];
    sigma[ip] = sigmaRho[ip]/rho[ip];
}

rho.updateSum();
//sigmaRho = rho.timesSigma();
sigmaRho.updateSum();
theta.updateSum();
sigmaRho.write();

Info << "Initialising sigmaRhof and sigmaf - values on the faces\n" << endl;
partitionedSurfaceScalarField sigmaRhof
(
    "sigmaRhof", partNames, fvc::interpolate(sigmaRho[0], "sigmaRho")
);
for(label ip = 1; ip < partNames.size(); ip++)
{
    sigmaRhof[ip] = fvc::interpolate(sigmaRho[ip], "sigmaRho");
}
sigmaRhof.updateSum();

partitionedSurfaceScalarField thetaf
(
    "thetaf", partNames, fvc::interpolate(theta[0]), sigmaRhof
);
for(label ip = 1; ip < partNames.size(); ip++)
{
    thetaf[ip] = fvc::interpolate(theta[ip]);
}
thetaf.updateSum();

partitionedSurfaceScalarField sigmaf
(
    "sigmaf", partNames, fvc::interpolate(sigma[0], "sigma")
);
for(label ip = 1; ip < partNames.size(); ip++)
{
    sigmaf[ip] = fvc::interpolate(sigma[ip], "sigma");
}
sigmaf.updateSum();

Info << "Reading in the velocity in each partition\n" << endl;
partitionedVolVectorField u
(
    "u", partNames, mesh, runTime.timeName(), sigmaRho
);
partitionedSurfaceVectorField Uf
(
    "Uf", partNames, fvc::interpolate(u[0], "Uf"), sigmaRhof
);
for(label ip = 0; ip < partNames.size(); ip++)
{
    Uf[ip] = fvc::interpolate(u[ip], "Uf");
}

/*partitionedSurfaceVectorField Uf
(
    "Uf", partNames, mesh, runTime.timeName(), sigmaRhof
);*/

// fixed flux boundary conditions
wordList fluxBCs(mesh.boundaryMesh().size(), "calculated");
forAll(mesh.boundaryMesh(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    if (btype == "wall" || btype == "symmetryPlane" || btype == "empty")
    {
        fluxBCs[patchi] = "fixedValue";
    }
}

Info << "Calculating the mass flux and volume flux in each partition\n"<<endl;
partitionedSurfaceScalarField volFlux
(
    "volFlux",
    partNames,
    surfaceScalarField
    (
        IOobject("volFlux"+'.'+partNames[0], runTime.timeName(), mesh),
        Uf[0] & mesh.Sf(),
        fluxBCs
    ),
    sigmaRhof
);
for(label ip = 1; ip < volFlux.size(); ip++)
{
    volFlux[ip] = Uf[ip] & mesh.Sf();
}
volFlux.updateSum();

partitionedSurfaceScalarField flux
(
    "flux", partNames, linearInterpolate(rho[0])*volFlux[0], sigmaf
);
for(label ip = 1; ip < volFlux.size(); ip++)
{
    flux[ip] = linearInterpolate(rho[ip])*volFlux[ip];
}
flux.updateSum();

partitionedSurfaceScalarField sigmaFlux = volFlux.timesSigma();

Info << "Reading in BCs and calculating the cell centre velocity, u, in each partition\n" << endl;
/*partitionedVolVectorField u
(
    "u", partNames, mesh, runTime.timeName(), sigmaRho
);
for(label ip = 0; ip < partNames.size(); ip++)
{
    u[ip] = fvc::reconstruct(volFlux[ip]);
}*/

partitionedVolScalarField w
(
    "w", partNames, u[0].component(vector::Z)
);
for(label ip = 0; ip < partNames.size(); ip++)
{
    w[ip] = u[ip].component(vector::Z);
}

partitionedVolScalarField wVar
(
    "wVar", partNames, mesh, runTime.timeName(), sigmaRho
);

Info << "Storing necessary old times and rates of change\n" << endl;
volFlux.storeTime();
theta.storeTime();
thetaVar.storeTime();
w.storeTime();
wVar.storeTime();
sigmaRho.storeTime();
flux.sum().oldTime();
flux.storeTime();
rho.storeTime();
sigma.storeTime();
sigmaFlux.storeTime();

partitionedVolScalarField sigmaRhoPreTransfer = sigmaRho;
partitionedVolScalarField thetaPreTransfer = theta;
partitionedVolScalarField thetaVarPreTransfer = thetaVar;
partitionedVolScalarField wPreTransfer = w;
partitionedVolScalarField wVarPreTransfer = wVar;

partitionedVolScalarField thetaPreTransferTot = theta;
partitionedVolScalarField thetaVarPreTransferTot = thetaVar;
partitionedVolScalarField wPreTransferTot = w;
partitionedVolScalarField wVarPreTransferTot = wVar;

Info << "The compressibility\n" << endl;
/*partitionedVolScalarField Psi
(
    "Psi",
    partNames,
    (1-kappa)/kappa*rho[0]/Exner,
    sigma
);
for (label ip = 1; ip < nParts; ip++)
{
    (1-kappa)/kappa*rho[ip]/Exner;
}
Psi.updateSum();*/
volScalarField Psi = (1-kappa)/kappa*rho.sum()/Exner;

// Alternative version of theta.sum
volScalarField rhoTheta = pRef/R*pow(Exner, (1-kappa)/kappa);

surfaceScalarField gradPcoeffSum
(
    IOobject("gradPcoeff", runTime.timeName(), mesh),
    Cp*fvc::interpolate(rhoTheta, "theta")
);

Info << "Initialising the drag between partitions" << endl;
surfaceScalarField cloudRadius
(
    "cloudRadius",
    max(cloudRadiusMin, cloudRadiusMax*sigmaf[0]*sigmaf[1])
);
surfaceScalarField dragCommon = -dragCoeff*sigmaRhof.sum()/cloudRadius
                         * mag(Uf[0] - Uf[1]);

#include "calculateDrag.H"

Info << "Initialising the mass transfer between partitions\n" << endl;
// Temporary mass and heat transfer fields (for calculations)
const dimensionedScalar massTransferZero
(
    "massTransferZero", dimensionSet(1,-3,-1,0,0), scalar(0)
);
const dimensionedScalar thetaTransferZero
(
    "massTransferZero", dimensionSet(0,0,-1,1,0), scalar(0)
);
const dimensionedScalar thetaVarTransferZero
(
    "massTransferZero", dimensionSet(0,0,-1,2,0), scalar(0)
);
const dimensionedScalar wTransferZero
(
    "wTransferZero", dimensionSet(0,1,-2,0,0), scalar(0)
);
const dimensionedScalar wVarTransferZero
(
    "wVarTransferZero", dimensionSet(0,2,-3,0,0), scalar(0)
);
const dimensionedScalar rateZero
(
    "rateZero", dimensionSet(0,0,-1,0,0), scalar(0)
);
const dimensionedScalar dimlessZero
(
    "dimlessZero", dimless, scalar(0)
);
const dimensionedScalar heatTransferZero
(
    "heatTransferZero", dimensionSet(0,-2,0,1,0), scalar(0)
);
const dimensionedScalar densityZero
(
    "densityZero", dimensionSet(1,-3,0,0,0), scalar(0)
);
const dimensionedScalar zeroTemp
(
    "zeroTemp", dimensionSet(0,0,0,1,0), scalar(0)
);
volScalarField massTransferTmp
(
    IOobject("massTransfer", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    massTransferZero
);
volScalarField thetaTransferTmp
(
    IOobject("thetaTransfer", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    thetaTransferZero
);
volScalarField thetaVarTransferTmp
(
    IOobject("thetaVarTransfer", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    thetaVarTransferZero
);
volScalarField wTransferTmp
(
    IOobject("wTransfer", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    wTransferZero
);
volScalarField wVarTransferTmp
(
    IOobject("wVarTransfer", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    wVarTransferZero
);

PtrList<partitionedVolScalarField> massTransfer(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    massTransfer.set
    (
        jp, 
        new partitionedVolScalarField
        (
            "massTransfer."+partNames[jp], partNames, massTransferTmp
        )
    );
}

PtrList<partitionedVolScalarField> massTransferAnalytic(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    massTransferAnalytic.set
    (
        jp, 
        new partitionedVolScalarField
        (
            "massTransferAnalytic."+partNames[jp], partNames, massTransferTmp
        )
    );
}

PtrList<partitionedVolScalarField> massTransferEntrainment(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    massTransferEntrainment.set
    (
        jp, 
        new partitionedVolScalarField
        (
            "massTransferEntrainment."+partNames[jp], partNames, massTransferTmp
        )
    );
}

PtrList<partitionedVolScalarField> massTransferDivTransfer(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    massTransferDivTransfer.set
    (
        jp, 
        new partitionedVolScalarField
        (
            "massTransferDivTransfer."+partNames[jp], partNames, massTransferTmp
        )
    );
}

PtrList<partitionedVolScalarField> massTransferSigmaDiffusion(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    massTransferSigmaDiffusion.set
    (
        jp, 
        new partitionedVolScalarField
        (
            "massTransferSigmaDiffusion."+partNames[jp], partNames, massTransferTmp
        )
    );
}

PtrList<partitionedVolScalarField> massTransferThetaDiffusion(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    massTransferThetaDiffusion.set
    (
        jp, 
        new partitionedVolScalarField
        (
            "massTransferThetaDiffusion."+partNames[jp], partNames, massTransferTmp
        )
    );
}

PtrList<partitionedVolScalarField> massTransferThetaVar(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    massTransferThetaVar.set
    (
        jp, 
        new partitionedVolScalarField
        (
            "massTransferThetaVar."+partNames[jp], partNames, massTransferTmp
        )
    );
}

PtrList<partitionedVolScalarField> massTransferWvar(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    massTransferWvar.set
    (
        jp, 
        new partitionedVolScalarField
        (
            "massTransferWvar."+partNames[jp], partNames, massTransferTmp
        )
    );
}


partitionedVolScalarField sigmaRhoAdvection
(
    "sigmaRhoAdvection", partNames, 0*sigmaRho[0]
);

partitionedVolScalarField thetaAdvection
(
    "thetaAdvection", partNames, 0*theta[0]
);

partitionedVolVectorField uAdvection
(
    "uAdvection", partNames, u[0]
);

partitionedVolVectorField uPressureGradient
(
    "uPressureGradient", partNames, u[0]
);

partitionedVolVectorField uPressureGradientPrime
(
    "uPressureGradientPrime", partNames, u[0]
);


partitionedVolScalarField thetaTransferChange
(
    "thetaTransferChange", partNames, 0*theta[0]/dt
);

partitionedVolScalarField thetaTransferAnalyticMean
(
    "thetaTransferAnalyticMean", partNames, thetaTransferChange[0]
);
partitionedVolScalarField thetaTransferAnalyticZero
(
    "thetaTransferAnalyticZero", partNames, thetaTransferChange[0]
);
partitionedVolScalarField thetaTransferAnalyticVar
(
    "thetaTransferAnalyticVar", partNames, thetaTransferChange[0]
);
partitionedVolScalarField thetaTransferAnalyticVarSmooth
(
    "thetaTransferAnalyticVarSmooth", partNames, thetaTransferChange[0]
);
partitionedVolScalarField thetaTransferAnalyticMeanVar
(
    "thetaTransferAnalyticMeanVar", partNames, thetaTransferChange[0]
);

partitionedVolScalarField thetaTransferDivMean
(
    "thetaTransferDivMean", partNames, thetaTransferChange[0]
);
partitionedVolScalarField thetaTransferDivVar
(
    "thetaTransferDivVar", partNames, thetaTransferChange[0]
);
partitionedVolScalarField thetaTransferDivMeanVar
(
    "thetaTransferDivMeanVar", partNames, thetaTransferChange[0]
);

partitionedVolScalarField thetaTransferSigmaDiffMean
(
    "thetaTransferSigmaDiffMean", partNames, thetaTransferChange[0]
);
partitionedVolScalarField thetaTransferSigmaDiffVar
(
    "thetaTransferSigmaDiffVar", partNames, thetaTransferChange[0]
);
partitionedVolScalarField thetaTransferSigmaDiffMeanVar
(
    "thetaTransferSigmaDiffMeanVar", partNames, thetaTransferChange[0]
);

partitionedVolScalarField thetaTransferThetaDiffMean
(
    "thetaTransferThetaDiffMean", partNames, thetaTransferChange[0]
);
partitionedVolScalarField thetaTransferThetaDiffVar
(
    "thetaTransferThetaDiffVar", partNames, thetaTransferChange[0]
);
partitionedVolScalarField thetaTransferThetaDiffMeanVar
(
    "thetaTransferThetaDiffMeanVar", partNames, thetaTransferChange[0]
);

partitionedVolScalarField thetaVarTransferChange
(
    "thetaVarTransferChange", partNames, 0*thetaVar[0]/dt
);

partitionedVolScalarField thetaVarTransferAnalyticMean
(
    "thetaVarTransferAnalyticMean", partNames, thetaVarTransferChange[0]
);
partitionedVolScalarField thetaVarTransferAnalyticZero
(
    "thetaVarTransferAnalyticZero", partNames, thetaVarTransferChange[0]
);
partitionedVolScalarField thetaVarTransferAnalyticVar
(
    "thetaVarTransferAnalyticVar", partNames, thetaVarTransferChange[0]
);
partitionedVolScalarField thetaVarTransferAnalyticVarSmooth
(
    "thetaVarTransferAnalyticVarSmooth", partNames, thetaVarTransferChange[0]
);
partitionedVolScalarField thetaVarTransferAnalyticMeanVar
(
    "thetaVarTransferAnalyticMeanVar", partNames, thetaVarTransferChange[0]
);

partitionedVolScalarField thetaVarTransferDivMean
(
    "thetaVarTransferDivMean", partNames, thetaVarTransferChange[0]
);
partitionedVolScalarField thetaVarTransferDivVar
(
    "thetaVarTransferDivVar", partNames, thetaVarTransferChange[0]
);
partitionedVolScalarField thetaVarTransferDivMeanVar
(
    "thetaVarTransferDivMeanVar", partNames, thetaVarTransferChange[0]
);

partitionedVolScalarField thetaVarTransferSigmaDiffMean
(
    "thetaVarTransferSigmaDiffMean", partNames, thetaVarTransferChange[0]
);
partitionedVolScalarField thetaVarTransferSigmaDiffVar
(
    "thetaVarTransferSigmaDiffVar", partNames, thetaVarTransferChange[0]
);
partitionedVolScalarField thetaVarTransferSigmaDiffMeanVar
(
    "thetaVarTransferSigmaDiffMeanVar", partNames, thetaVarTransferChange[0]
);

partitionedVolScalarField thetaVarTransferThetaDiffMean
(
    "thetaVarTransferThetaDiffMean", partNames, thetaVarTransferChange[0]
);
partitionedVolScalarField thetaVarTransferThetaDiffVar
(
    "thetaVarTransferThetaDiffVar", partNames, thetaVarTransferChange[0]
);
partitionedVolScalarField thetaVarTransferThetaDiffMeanVar
(
    "thetaVarTransferThetaDiffMeanVar", partNames, thetaVarTransferChange[0]
);


partitionedVolScalarField wTransferChange
(
    "wTransferChange", partNames, 0*w[0]/dt
);

partitionedVolScalarField wTransferAnalyticMean
(
    "wTransferAnalyticMean", partNames, wTransferChange[0]
);
partitionedVolScalarField wTransferAnalyticVar
(
    "wTransferAnalyticVar", partNames, wTransferChange[0]
);
partitionedVolScalarField wTransferAnalyticVarSmooth
(
    "wTransferAnalyticVarSmooth", partNames, wTransferChange[0]
);
partitionedVolScalarField wTransferAnalyticMeanVar
(
    "wTransferAnalyticMeanVar", partNames, wTransferChange[0]
);
partitionedVolScalarField wTransferAnalyticZero
(
    "wTransferAnalyticZero", partNames, wTransferChange[0]
);

partitionedVolScalarField wTransferDivMean
(
    "wTransferDivMean", partNames, wTransferChange[0]
);
partitionedVolScalarField wTransferDivVar
(
    "wTransferDivVar", partNames, wTransferChange[0]
);
partitionedVolScalarField wTransferDivMeanVar
(
    "wTransferDivMeanVar", partNames, wTransferChange[0]
);
partitionedVolScalarField wTransferDivZero
(
    "wTransferDivZero", partNames, wTransferChange[0]
);

partitionedVolScalarField wTransferSigmaDiffMean
(
    "wTransferSigmaDiffMean", partNames, wTransferChange[0]
);
partitionedVolScalarField wTransferSigmaDiffVar
(
    "wTransferSigmaDiffVar", partNames, wTransferChange[0]
);
partitionedVolScalarField wTransferSigmaDiffMeanVar
(
    "wTransferSigmaDiffMeanVar", partNames, wTransferChange[0]
);
partitionedVolScalarField wTransferSigmaDiffZero
(
    "wTransferSigmaDiffZero", partNames, wTransferChange[0]
);

partitionedVolScalarField wTransferThetaDiffMean
(
    "wTransferThetaDiffMean", partNames, wTransferChange[0]
);
partitionedVolScalarField wTransferThetaDiffVar
(
    "wTransferThetaDiffVar", partNames, wTransferChange[0]
);
partitionedVolScalarField wTransferThetaDiffMeanVar
(
    "wTransferThetaDiffMeanVar", partNames, wTransferChange[0]
);
partitionedVolScalarField wTransferThetaDiffZero
(
    "wTransferThetaDiffZero", partNames, wTransferChange[0]
);


partitionedVolScalarField wVarTransferChange
(
    "wVarTransferChange", partNames, 0*wVar[0]/dt
);

partitionedVolScalarField wVarTransferAnalyticMean
(
    "wVarTransferAnalyticMean", partNames, wVarTransferChange[0]
);
partitionedVolScalarField wVarTransferAnalyticVar
(
    "wVarTransferAnalyticVar", partNames, wVarTransferChange[0]
);
partitionedVolScalarField wVarTransferAnalyticVarSmooth
(
    "wVarTransferAnalyticVarSmooth", partNames, wVarTransferChange[0]
);
partitionedVolScalarField wVarTransferAnalyticMeanVar
(
    "wVarTransferAnalyticMeanVar", partNames, wVarTransferChange[0]
);
partitionedVolScalarField wVarTransferAnalyticZero
(
    "wVarTransferAnalyticZero", partNames, wVarTransferChange[0]
);

partitionedVolScalarField wVarTransferDivMean
(
    "wVarTransferDivMean", partNames, wVarTransferChange[0]
);
partitionedVolScalarField wVarTransferDivVar
(
    "wVarTransferDivVar", partNames, wVarTransferChange[0]
);
partitionedVolScalarField wVarTransferDivMeanVar
(
    "wVarTransferDivMeanVar", partNames, wVarTransferChange[0]
);
partitionedVolScalarField wVarTransferDivZero
(
    "wVarTransferDivZero", partNames, wVarTransferChange[0]
);

partitionedVolScalarField wVarTransferSigmaDiffMean
(
    "wVarTransferSigmaDiffMean", partNames, wVarTransferChange[0]
);
partitionedVolScalarField wVarTransferSigmaDiffVar
(
    "wVarTransferSigmaDiffVar", partNames, wVarTransferChange[0]
);
partitionedVolScalarField wVarTransferSigmaDiffMeanVar
(
    "wVarTransferSigmaDiffMeanVar", partNames, wVarTransferChange[0]
);
partitionedVolScalarField wVarTransferSigmaDiffZero
(
    "wVarTransferSigmaDiffZero", partNames, wVarTransferChange[0]
);

partitionedVolScalarField wVarTransferThetaDiffMean
(
    "wVarTransferThetaDiffMean", partNames, wVarTransferChange[0]
);
partitionedVolScalarField wVarTransferThetaDiffVar
(
    "wVarTransferThetaDiffVar", partNames, wVarTransferChange[0]
);
partitionedVolScalarField wVarTransferThetaDiffMeanVar
(
    "wVarTransferThetaDiffMeanVar", partNames, wVarTransferChange[0]
);
partitionedVolScalarField wVarTransferThetaDiffZero
(
    "wVarTransferThetaDiffZero", partNames, wVarTransferChange[0]
);


partitionedVolScalarField wTransferGradDiv
(
    "wTransferGradDiv", partNames, 0*w[0]/dt
);

partitionedVolScalarField wTransferDrag
(
    "wTransferDrag", partNames, 0*w[0]/dt
);


partitionedVolScalarField thetaT
(
    "thetaT", partNames, 0*theta[0]
);
partitionedVolScalarField thetaTAnalyticMean
(
    "thetaTAnalyticMean", partNames, 0*theta[0]
);
partitionedVolScalarField thetaTAnalyticZero
(
    "thetaTAnalyticZero", partNames, 0*theta[0]
);
partitionedVolScalarField thetaTAnalyticVar
(
    "thetaTAnalyticVar", partNames, 0*theta[0]
);
partitionedVolScalarField thetaTAnalyticVarSmooth
(
    "thetaTAnalyticVarSmooth", partNames, 0*theta[0]
);
partitionedVolScalarField thetaTAnalyticMeanVar
(
    "thetaTAnalyticMeanVar", partNames, 0*theta[0]
);
partitionedVolScalarField thetaVarT
(
    "thetaVarT", partNames, 0*thetaVar[0]
);

partitionedVolScalarField wT
(
    "wT", partNames, 0*w[0]
);
partitionedVolScalarField wTAnalyticMean
(
    "wTAnalyticMean", partNames, 0*w[0]
);
partitionedVolScalarField wTAnalyticZero
(
    "wTAnalyticZero", partNames, 0*w[0]
);
partitionedVolScalarField wTAnalyticVar
(
    "wTAnalyticVar", partNames, 0*w[0]
);
partitionedVolScalarField wTAnalyticVarSmooth
(
    "wTAnalyticVarSmooth", partNames, 0*w[0]
);
partitionedVolScalarField wTAnalyticMeanVar
(
    "wTAnalyticMeanVar", partNames, 0*w[0]
);
partitionedVolScalarField wVarT
(
    "wVarT", partNames, 0*wVar[0]
);


Info << "Initialising the theta transfer between partitions\n" << endl;
const dimensionedScalar thetaTransferTermZero
(
    "thetaTransferTermZero", dimensionSet(0,0,-1,0,0), scalar(0)
);
volScalarField thetaTransferTermTmp
(
    IOobject("thetaTransferTerm", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    thetaTransferTermZero
);
PtrList<partitionedVolScalarField> thetaTransferTerm(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    thetaTransferTerm.set
    (
        jp, 
        new partitionedVolScalarField
        (
            partNames[jp]+'.'+"thetaTransferTerm", partNames, thetaTransferTermTmp
        )
    );
}

const dimensionedScalar thetaVarTransferTermZero
(
    "thetaVarTransferTermZero", dimensionSet(0,0,-1,-1,0), scalar(0)
);
volScalarField thetaVarTransferTermTmp
(
    IOobject("thetaVarTransferTerm", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    thetaVarTransferTermZero
);
PtrList<partitionedVolScalarField> thetaVarTransferTerm(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    thetaVarTransferTerm.set
    (
        jp, 
        new partitionedVolScalarField
        (
            partNames[jp]+'.'+"thetaVarTransferTerm", partNames, thetaVarTransferTermTmp
        )
    );
}

Info << "Initialising the w transfer between partitions\n" << endl;
const dimensionedScalar wTransferTermZero
(
    "wTransferTermZero", dimensionSet(0,0,-1,0,0), scalar(0)
);
volScalarField wTransferTermTmp
(
    IOobject("wTransferTerm", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    wTransferTermZero
);
PtrList<partitionedVolScalarField> wTransferTerm(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    wTransferTerm.set
    (
        jp, 
        new partitionedVolScalarField
        (
            partNames[jp]+'.'+"wTransferTerm", partNames, wTransferTermTmp
        )
    );
}

const dimensionedScalar wVarTransferTermZero
(
    "wVarTransferTermZero", dimensionSet(0,-1,0,0,0), scalar(0)
);
volScalarField wVarTransferTermTmp
(
    IOobject("wVarTransferTerm", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    wVarTransferTermZero
);
PtrList<partitionedVolScalarField> wVarTransferTerm(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    wVarTransferTerm.set
    (
        jp, 
        new partitionedVolScalarField
        (
            partNames[jp]+'.'+"wVarTransferTerm", partNames, wVarTransferTermTmp
        )
    );
}

// Transfer terms for other variables (assuming nParts == 2)
partitionedVolScalarField transfer
(
    "transfer",
    partNames,
    volScalarField
    (
        IOobject("transfer", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("", dimless, scalar(0))
    )
);

partitionedSurfaceScalarField transferf
(
    "transferf",
    partNames,
    surfaceScalarField
    (
        IOobject("transfer", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("", dimless, scalar(0))
    )
);

partitionedVolScalarField transferRate
(
    "transferRate",
    partNames,
    volScalarField
    (
        IOobject("transferRate", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("", dimensionSet(0,0,-1,0,0), scalar(0))
    )
);

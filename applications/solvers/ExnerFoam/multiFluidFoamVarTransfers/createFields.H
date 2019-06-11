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

Info << "Setting rho for each partition\n" << endl;
partitionedVolScalarField rho
(
    "rho",
    partNames,
//    mesh, runTime.timeName(),
    volScalarField
    (
        IOobject("rho", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("rho", dimDensity, scalar(1)),
        sigma[0].boundaryField().types()
    ),
    sigma
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
for(label ip = 0; ip < partNames.size(); ip++)
{
    rho[ip] = pRef/(R*theta[ip])*pow(Exner, (1-kappa)/kappa);
}
rho.updateSum();
sigmaRho = rho.timesSigma();
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
partitionedSurfaceVectorField Uf
(
    "Uf", partNames, mesh, runTime.timeName(), sigmaRhof
);

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
partitionedVolVectorField u
(
    "u", partNames, mesh, runTime.timeName(), sigmaRho
);
for(label ip = 0; ip < partNames.size(); ip++)
{
    u[ip] = fvc::reconstruct(volFlux[ip]);
}

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
    "wVar", partNames, 0*sqr(w[0])
);

partitionedVolScalarField KwVar
(
    "KwVar", partNames, dimensionedScalar("", dimensionSet(0,2,-1,0,0), scalar(0))*sigma[0]
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


//Additional coefficients
const dimensionedScalar sigmaRhoMin("",sigmaRho[0].dimensions(),scalar(1e-16));
const dimensionedScalar thetaMin("",theta[0].dimensions(),scalar(1e-16));
const dimensionedScalar thetaVarMin("",thetaVar[0].dimensions(),scalar(1e-16));
const dimensionedScalar wMin("",w[0].dimensions(),scalar(1e-16));
const dimensionedScalar wVarMin("",wVar[0].dimensions(),scalar(1e-16));
const dimensionedScalar wOne("",w[0].dimensions(),scalar(1));

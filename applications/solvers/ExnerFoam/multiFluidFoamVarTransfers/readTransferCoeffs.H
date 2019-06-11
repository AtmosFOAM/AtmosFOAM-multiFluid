Info << "\nReading environmentalProperties" << endl;

IOdictionary transferProps
(
    IOobject
    (
        "transferProperties",
        runTime.system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const dimensionedScalar cloudRadiusMax(transferProps.lookup("cloudRadiusMax"));
const dimensionedScalar cloudRadiusMin(transferProps.lookup("cloudRadiusMin"));
const scalar dragCoeff(readScalar(transferProps.lookup("dragCoeff")));
const scalar divCoeff(readScalar(transferProps.lookup("divCoeff")));

// Diffusion coefficient for inter-partition mass transfer
const dimensionedScalar Ksigma(transferProps.lookup("Ksigma"));
const dimensionedScalar Kw(transferProps.lookup("Kw"));
const dimensionedScalar Ktheta(transferProps.lookup("Ktheta"));

// Is mass transferred based on partition divergence
const Switch divTransfer(transferProps.lookup("divTransfer"));

// Is mass transferred based on mixing (drag) at interfaces
const Switch dragTransfer(transferProps.lookup("dragTransfer"));
const scalar dragTransferCoeff(readScalar(transferProps.lookup("dragTransferCoeff")));

// Is mass transfer decided based on vertical velocity
const Switch wTransfer(transferProps.lookup("wTransfer"));

// Is mass transfer decided based on grad(theta)
const Switch thetaTransfer(transferProps.lookup("thetaTransfer"));
const dimensionedScalar thetaTransferDiffusivity
(
    transferProps.lookup("thetaTransferDiffusivity")
);

// Is heat transferred based on overlapping theta profiles
const Switch thetaVarTransfer(transferProps.lookup("thetaVarTransfer"));
const dimensionedScalar thetaVarTimescale
(
    transferProps.lookup("thetaVarTimescale")
);

const Switch wVarTransfer(transferProps.lookup("wVarTransfer"));
const dimensionedScalar wVarTimescale
(
    transferProps.lookup("wVarTimescale")
);


// Only do heat transfers where fluids have been mixed?
const Switch localThetaVarTransfer(transferProps.lookup("localThetaVarTransfer"));
const Switch varMassTransfer(transferProps.lookup("varMassTransfer"));

// Minimum sigmaRho for each partition
const scalar minSigma = readScalar(transferProps.lookup("minSigma"));
const scalar minSigmaDiv = max(minSigma, VSMALL);

IOdictionary blockMeshDict
(
    IOobject
    (
        "blockMeshDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const scalar xmin = readScalar(blockMeshDict.lookup("xmin"));
const scalar xmax = readScalar(blockMeshDict.lookup("xmax"));
const scalar nx = readScalar(blockMeshDict.lookup("nx"));
const scalar convertToMeters = readScalar(blockMeshDict.lookup("convertToMeters"));
const scalar deltax = (xmax-xmin)/nx * convertToMeters;
scalar thetaVarLengthScale = min(1, deltax/10000.);

// Diffusion based on velocity variance (accounts for sub-filter min/max speeds)
const Switch wVarDiffusion(transferProps.lookup("wVarDiffusion"));
const dimensionedScalar KwVariance(transferProps.lookup("KwVariance"));
const dimensionedScalar dz("dz", dimLength, scalar(100));




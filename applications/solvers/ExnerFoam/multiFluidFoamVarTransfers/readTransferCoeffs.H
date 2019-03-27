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

// Diffusion coefficient for inter-partition mass transfer
const dimensionedScalar Ksigma(transferProps.lookup("Ksigma"));
const dimensionedScalar Kw(transferProps.lookup("Kw"));
const dimensionedScalar Ktheta(transferProps.lookup("Ktheta"));

// Is mass transferred based on partition divergence
const Switch divTransfer(transferProps.lookup("divTransfer"));

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
const Switch wVarTransfer(transferProps.lookup("wVarTransfer"));
const dimensionedScalar thetaVarTimescale
(
    transferProps.lookup("thetaVarTimescale")
);

// Change velocities proportionally to amount temperatures change.
const Switch thetaVarTransferWcoupling(transferProps.lookup("thetaVarTransferWcoupling"));
const scalar wCouplingCoeff = readScalar(transferProps.lookup("wCouplingCoeff"));

// Only do heat transfers where fluids have been mixed?
const Switch localThetaVarTransfer(transferProps.lookup("localThetaVarTransfer"));

// Minimum sigmaRho for each partition
const scalar minSigma = readScalar(transferProps.lookup("minSigma"));
const scalar minSigmaDiv = max(minSigma, VSMALL);
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

// Is mass transferred based on partition divergence
const Switch divTransfer(transferProps.lookup("divTransfer"));

// Is mass transfer decided based on vertical velocity
const Switch wTransfer(transferProps.lookup("wTransfer"));
const dimensionedScalar transferTimeScale = wTransfer ? 
      dimensionedScalar(transferProps.lookup("transferTimeScale")) :
      dimensionedScalar("transferProps", dimTime, scalar(0));

// Is mass transfer decided based on grad(theta)
const Switch thetaTransfer(transferProps.lookup("thetaTransfer"));
const dimensionedScalar thetaTransferDiffusivity
(
    transferProps.lookup("thetaTransferDiffusivity")
);
const dimensionedScalar thetaTransferDiffusivity2
(
    transferProps.lookupOrDefault<dimensionedScalar>
    (
        "thetaTransferDiffusivity2",
        dimensionedScalar("zero", thetaTransferDiffusivity.dimensions(), scalar(0))
    )
);

// Minimum sigmaRho for each partition
const scalar minSigma = readScalar(transferProps.lookup("minSigma"));

// Max and min temperature perturbation - the difference in temperature of
// the transferred variable. thetat[0] > 0 and thetat[1] < 0
const dimensionedScalar thetatMin
(
    "thetatMin",
    dimTemperature,
    transferProps.lookupOrDefault<scalar>("thetatMin", 0)
);
const dimensionedScalar thetatMax
(
    "thetatMax",
    dimTemperature,
    transferProps.lookupOrDefault<scalar>("thetatMax", 0)
);


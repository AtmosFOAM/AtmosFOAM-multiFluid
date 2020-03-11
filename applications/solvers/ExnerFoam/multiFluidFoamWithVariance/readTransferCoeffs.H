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

// Minimum sigma for each partition
const scalar minSigma = readScalar(transferProps.lookup("minSigma"));

// Transfer to keep flagged cells uniform sigma
const Switch uniformSigma(transferProps.lookup("uniformSigma"));

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


// Parameters used for variance variables
// Use temperature and/or velocity variance
const Switch useThetaVariance
(
    transferProps.lookupOrDefault<bool>("useThetaVariance", false)
);

// Does the transferred temperature profile 
// have a Smooth profile or Sharp boundary?
const word thetaVarTransferType
(
    transferProps.lookupOrDefault<word>("thetaVarTransferType", "Sharp")
);

// Minimum value for divisions by variance
const dimensionedScalar thetaVarMin("",dimensionSet(0,0,0,2,0),scalar(1e-16));

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

// Diffusion coefficient for inter-partition mass transfer
const dimensionedScalar Ksigma(transferProps.lookup("Ksigma"));

// Is mass transferred based on partition divergence
const Switch divTransfer(transferProps.lookup("divTransfer"));

// Is mass transfer decided based on vertical velocity
const Switch wTransfer(transferProps.lookup("wTransfer"));
const dimensionedScalar transferTimeScale = wTransfer ? 
      dimensionedScalar(transferProps.lookup("transferTimeScale")) :
      dimensionedScalar("transferProps", dimTime, scalar(0));

// Is mass transfer decided based on grad(b)
const Switch buoyancyTransfer(transferProps.lookup("buoyancyTransfer"));
const dimensionedScalar bTransferDiffusivity
(
    transferProps.lookup("bTransferDiffusivity")
);

// Minimum sigma for each partition
const scalar minSigma = readScalar(transferProps.lookup("minSigma"));

// Transfer to keep flagged cells uniform sigma
const Switch uniformSigma(transferProps.lookup("uniformSigma"));


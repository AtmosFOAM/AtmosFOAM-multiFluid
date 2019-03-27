Info << "\nReading transferProperties" << endl;

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

// The amount of divergence transfer used (between zero and one)
const scalar divTransfer(readScalar(transferProps.lookup("divTransfer")));
if (divTransfer < 0 || divTransfer > 1)
{
    FatalErrorIn("multiFluidBoussinesqFoam")
        << "divTransfer read in from file system/transferProperties should be "
        << "between zero and one, not " << divTransfer << exit(FatalError);
}

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


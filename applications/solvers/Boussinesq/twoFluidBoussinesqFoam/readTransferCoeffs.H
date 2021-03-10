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

// Drag parameters
const dimensionedScalar cloudRadius(transferProps.lookup("cloudRadius"));
const scalar dragCoeff(readScalar(transferProps.lookup("dragCoeff")));

// Note if there are no transfers from any sources
const Switch noTransfers(transferProps.lookup("noTransfers"));

// Diffusion coefficient for inter-partition mass transfer
const dimensionedScalar Ksigma(transferProps.lookup("Ksigma"));

// The amount of divergence transfer used (between zero and one)
const scalar divTransfer = noTransfers? scalar(0) :
      readScalar(transferProps.lookup("divTransfer"));
if (divTransfer < 0 || divTransfer > 2)
{
    FatalErrorIn("multiFluidBoussinesqFoam")
        << "divTransfer read in from file system/transferProperties should be "
        << "between zero and one, not " << divTransfer << exit(FatalError);
}

// Transfer based on dwdz
const scalar dwdzTransfer = noTransfers ? scalar(0) :
    readScalar(transferProps.lookup("dwdzTransfer"));

// Transfers based on w
const scalar wTransfer = readScalar(transferProps.lookup("wTransfer"));

// Transfer w=0 (w^T = 0) (scaled by sigma)
const Switch wTransfer0 = noTransfers? Switch(false) :
      transferProps.lookup("wTransfer0");

// Proportion of sub-grid-scale variability of w
const scalar wPrimeCoeff = noTransfers? scalar(0) :
      readScalar(transferProps.lookup("wPrimeCoeff"));

// Minimum sigma for each partition
const scalar minSigma = noTransfers? scalar(0) :
      readScalar(transferProps.lookup("minSigma"));

// Transfer to keep flagged cells uniform sigma
const Switch uniformSigma = noTransfers? Switch(false) :
      Switch(transferProps.lookup("uniformSigma"));

// Transfer buoyancy anomalies based on linear combination of mean fluid 
// buoyancies
const scalar bAnomalyLinearCoeff = transferProps.lookupOrDefault<scalar>
(
    "bAnomalyLinearCoeff", scalar(0)
);

// Transfer buoyancy anomalies based on constant factor times mean buoyancy
const scalar bAnomalyMeanCoeff = transferProps.lookupOrDefault<scalar>
(
    "bAnomalyMeanCoeff", scalar(0)
);

// Transfer buoyancy to keep fluid 0 decending and fluid 1 ascending
const Switch keepWseparate = transferProps.lookup("keepWseparate");

// Transfer based on laplacian of buoyancy
const dimensionedScalar bTransferDiffusivity
(
    transferProps.lookupOrDefault<dimensionedScalar>
    (
        "bTransferDiffusivity",
        dimensionedScalar("bTransferDiffusivity", dimensionSet(0,2,-1,0,0), scalar(0))
    )
);
const dimensionedScalar bTransferTimeScale
(
    transferProps.lookupOrDefault<dimensionedScalar>
    (
        "bTransferTimeScale",
        dimensionedScalar("bTransferTimeScale", dimTime, scalar(0))
    )
);
const dimensionedScalar bTransferMax
(
    transferProps.lookupOrDefault<dimensionedScalar>
    (
        "bTransferMax",
        dimensionedScalar("bTransferMax", dimensionSet(0,1,-2,0,0), scalar(0))
    )
);


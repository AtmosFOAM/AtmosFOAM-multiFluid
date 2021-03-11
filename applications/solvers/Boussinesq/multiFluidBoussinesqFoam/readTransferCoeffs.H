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

// The partition names
const wordList partNames(transferProps.lookup("fluids"));
const label nParts = partNames.size();

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

// Transfer based on dbdz
const scalar dbdzTransfer = noTransfers ? scalar(0) :
    readScalar(transferProps.lookup("dbdzTransfer"));

// Transfers based on w
const scalar wTransfer = readScalar(transferProps.lookup("wTransfer"));

// Transfer w=0 (w^T = 0) (scaled by sigma)
const Switch wTransfer0 = noTransfers? Switch(false) :
      transferProps.lookup("wTransfer0");

// Proportion of sub-grid-scale variability of w
const scalar wPrimeCoeff = noTransfers? scalar(0) :
      readScalar(transferProps.lookup("wPrimeCoeff"));

// Minimum sigma for each partition
const scalar minSigma = readScalar(transferProps.lookup("minSigma"));

// Transfer buoyancy anomalies based on constant factor times mean buoyancy
const scalar bAnomalyMeanCoeff = transferProps.lookupOrDefault<scalar>
(
    "bAnomalyMeanCoeff", scalar(0)
);


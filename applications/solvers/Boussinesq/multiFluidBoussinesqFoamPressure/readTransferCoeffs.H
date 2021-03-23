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

// Model for transfers between fluids
enum TransferType
{
    noTransfer, divTransfer, specialTransfer, diffusionOnlyTransfer, badTransfer
};
const word transferTypeWord(transferProps.lookup("transferType"));
const TransferType transferType = 
       transferTypeWord == "noTransfer" ? noTransfer :
       transferTypeWord == "divTransfer" ? divTransfer :
       transferTypeWord == "special" ? specialTransfer :
       transferTypeWord == "diffusionOnly" ? diffusionOnlyTransfer :
       badTransfer;

// Diffusion coefficient for inter-partition mass transfer
const dimensionedScalar Ksigma(transferProps.lookup("Ksigma"));

// Transfer w=0 (w^T = 0) (scaled by sigma)
const scalar wTransfer0 = transferType == noTransfer? scalar(0) :
      readScalar(transferProps.lookup("wTransfer0"));

// Minimum sigma for each partition and for partition 0
const scalarList minSigma(transferProps.lookup("minSigma"));

// Transfer buoyancy anomalies based on constant factor times mean buoyancy
const scalar bAnomalyMeanCoeff = readScalar
(
    transferProps.lookup("bAnomalyMeanCoeff")
);

// Check that Ksigma is non-zero if alpha or nu are non-zero
if
(
    (mag(alpha.value() - Ksigma.value()) > SMALL)
 && (mag(Ksigma.value()) < SMALL)
 && nParts > 1
)
{
    WarningIn("multiFluidBoussinesqFoam")
         << " if one out of alpha and Ksigma "
         << "is zero then they must both be zero but alpha = "
         << alpha.value() << " and Ksigma = " << Ksigma.value()
         << endl;
}
if
(
    (mag(nu.value() - Ksigma.value()) > SMALL)
 && (mag(Ksigma.value()) < SMALL)
 && nParts > 1
)
{
    WarningIn("multiFluidBoussinesqFoam")
        << " if one out of nu and Ksigma "
         << "is zero then they must both be zero but nu = "
         << nu.value() << " and Ksigma = " << Ksigma.value()
         << endl;
}


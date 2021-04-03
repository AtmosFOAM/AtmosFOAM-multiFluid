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

// Model for transfers between fluids
enum TransferType
{
    noTransfer, dwdzTransfer, divTransfer, specialTransfer, diffusionOnlyTransfer, badTransfer
};
const word transferTypeWord(transferProps.lookup("transferType"));
const TransferType transferType = 
       transferTypeWord == "noTransfer" ? noTransfer :
       transferTypeWord == "dwdzTransfer" ? dwdzTransfer:
       transferTypeWord == "divTransfer" ? divTransfer:
       transferTypeWord == "special" ? specialTransfer :
       transferTypeWord == "diffusionOnly" ? diffusionOnlyTransfer :
       badTransfer;

if (transferType == badTransfer)
{
    FatalErrorIn("readTransferCoeffs.H") 
             << " transferType must be one of "
             << "noTransfer, dwdzTransfer, divTransfer, specialTransfer, diffusionOnlyTransfer"
             << " not " << transferTypeWord << " in system/transferProperties"
             << exit(FatalError);
}

// Scale the mass transfes
const scalar transferScale = transferType == noTransfer? scalar(0) :
      readScalar(transferProps.lookup("transferScale"));

// Minimum sigma for each partition and for partition 0
const scalarList minSigma(transferProps.lookup("minSigma"));

// Transfer w=0 (w^T = 0) (scaled by sigma)
const scalar wTransferScale = transferType == noTransfer? scalar(0) :
      readScalar(transferProps.lookup("wTransferScale"));

// Transfer buoyancy anomalies based on constant factor times mean buoyancy
const scalar bTransferScale = readScalar
(
    transferProps.lookup("bTransferScale")
);

// Constant of proportionality between pressure and divergence
uniformDimensionedScalarField gamma
(
    IOobject("gamma", runTime.timeName(), mesh), 
    dimensionedScalar(transferProps.lookup("gamma"))
);

// Drag parameters
const dimensionedScalar cloudRadius(transferProps.lookup("cloudRadius"));
const scalar dragCoeff(readScalar(transferProps.lookup("dragCoeff")));


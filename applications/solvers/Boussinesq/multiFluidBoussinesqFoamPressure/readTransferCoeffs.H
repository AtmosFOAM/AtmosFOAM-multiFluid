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
    noTransfer, divTransfer, pressureTransfer, specialTransfer,
    diffusionOnlyTransfer, badTransfer
};
const word transferTypeWord(transferProps.lookup("transferType"));
const TransferType transferType = 
       transferTypeWord == "noTransfer" ? noTransfer :
       transferTypeWord == "divTransfer" ? divTransfer :
       transferTypeWord == "pressureTransfer" ? pressureTransfer:
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

// gamma for pressureTransfer
const dimensionedScalar gamma = transferType == pressureTransfer ?
      dimensionedScalar(transferProps.lookup("gamma")) :
      dimensionedScalar("gamma", dimensionSet(0,2,-1,0,0), scalar(0));


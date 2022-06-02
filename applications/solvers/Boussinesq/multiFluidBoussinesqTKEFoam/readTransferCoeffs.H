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
    noTransfer, dwdzTransfer, diffusionOnlyTransfer, badTransfer
};
const word transferTypeWord(transferProps.lookup("transferType"));
const TransferType transferType = 
       transferTypeWord == "noTransfer" ? noTransfer :
       transferTypeWord == "dwdzTransfer" ? dwdzTransfer:
       transferTypeWord == "diffusionOnly" ? diffusionOnlyTransfer :
       badTransfer;

if (transferType == badTransfer)
{
    FatalErrorIn("readTransferCoeffs.H") 
             << " transferType must be one of "
             << "noTransfer, dwdzTransfer, diffusionOnlyTransfer"
             << " not " << transferTypeWord << " in system/transferProperties"
             << exit(FatalError);
}

if (transferType == dwdzTransfer)
{
    if (nParts == 2 && (partNames[0] != "down" || partNames[1] != "up"))
    {
        FatalErrorIn("readTransferCoeffs.H") << " fluids " << partNames
             << " given but two fluids and transfer type dwdzTransfer the fluids must be down and up" << exit(FatalError);
    }
    else if
    (
        nParts == 3
     && (partNames[0]!="stable" || partNames[1]!="up" || partNames[2]!="down")
    )
    {
        FatalErrorIn("readTransferCoeffs.H") << " fluids " << partNames
             << " given but three fluids and transfer type dwdzTransfer the fluids must be stable, up and down"
             << exit(FatalError);
    }
}

// Scale the mixing transfer
const scalar mixingTransfer = transferType == noTransfer? scalar(0) :
                  readScalar(transferProps.lookup("mixingTransfer"));

// Scale the mass transfers
const scalar transferScale = transferType == noTransfer? scalar(0) :
      readScalar(transferProps.lookup("transferScale"));
      
// Entrainment and detrainment rates (both based on dwdz)
const scalar entRate = readScalar(transferProps.lookup("entRate"));
const scalar detrainRate = readScalar(transferProps.lookup("detrainRate"));

// Minimum sigma for each partition and for partition 0
const scalarList minSigma(transferProps.lookup("minSigma"));

// Vertical velocity transferred, w^T scaled by TKE
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


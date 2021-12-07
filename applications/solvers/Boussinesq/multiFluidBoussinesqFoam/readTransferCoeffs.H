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
    noTransfer, dwdzTransfer, dbdzTransfer, wTransfer, diffusionOnlyTransfer, badTransfer
};
const word transferTypeWord(transferProps.lookup("transferType"));
const TransferType transferType = 
       transferTypeWord == "noTransfer" ? noTransfer :
       transferTypeWord == "dwdzTransfer" ? dwdzTransfer:
       transferTypeWord == "dbdzTransfer" ? dbdzTransfer:
       transferTypeWord == "wTransfer" ? wTransfer :
       transferTypeWord == "diffusionOnly" ? diffusionOnlyTransfer :
       badTransfer;

if (transferType == badTransfer)
{
    FatalErrorIn("readTransferCoeffs.H") 
             << " transferType must be one of "
             << "noTransfer, dwdzTransfer, dbdzTransfer, wTransfer, diffusionOnlyTransfer"
             << " not " << transferTypeWord << " in system/transferProperties"
             << exit(FatalError);
}

if (transferType == dwdzTransfer || transferType == dbdzTransfer || transferType == wTransfer)
{
    if (nParts == 2 && (partNames[0] != "down" || partNames[1] != "up"))
    {
        FatalErrorIn("readTransferCoeffs.H") << " fluids " << partNames
             << " given but two fluids and transfer type wTransfer, dbdzTransfer or dwdzTransfer the fluids must be down and up" << exit(FatalError);
    }
    else if
    (
        nParts == 3
     && (partNames[0]!="stable" || partNames[1]!="up" || partNames[2]!="down")
    )
    {
        FatalErrorIn("readTransferCoeffs.H") << " fluids " << partNames
             << " given but three fluids and transfer type wTransfer, dbdzTransfer or dwdzTransfer the fluids must be stable, up and down"
             << exit(FatalError);
    }
}

// Smooth the mass transfers
const scalar massTransferSmooth
     = transferProps.lookupOrDefault<scalar>("massTransferSmooth", scalar(0));

// Scale the mass transfers
const scalar transferScale = transferType == noTransfer? scalar(0) :
      readScalar(transferProps.lookup("transferScale"));
      
// Entrainment and detrainment rates (both based on dwdz)
const scalar entRate = readScalar(transferProps.lookup("entRate"));
const scalar detrainRate = readScalar(transferProps.lookup("detrainRate"));

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

// Diffusion coefficient for inter-partition mass transfer
const dimensionedScalar Ksigma(transferProps.lookup("Ksigma"));

// Constant of proportionality between pressure and divergence
uniformDimensionedScalarField gamma
(
    IOobject("gamma", runTime.timeName(), mesh), 
    dimensionedScalar(transferProps.lookup("gamma"))
);

// Drag parameters
const dimensionedScalar cloudRadius(transferProps.lookup("cloudRadius"));
const scalar dragCoeff(readScalar(transferProps.lookup("dragCoeff")));

// Additional thermal diffusivity for fluid 0
// Scaling to go from distance to wall to thermal diffusivity of fluid 0
const dimensionedScalar alphaScale
(
    transferProps.lookupOrDefault<dimensionedScalar>
    (
        "alphaScale", 
        dimensionedScalar("alphaScale", dimVelocity, scalar(0))
    )
);
// Maximum thermal diffusivity of fluid 0
const dimensionedScalar maxAlpha0
(
    transferProps.lookupOrDefault<dimensionedScalar>
    (
        "maxAlpha0",
        dimensionedScalar("maxAlpha0", alpha.dimensions(), scalar(0))
    )
);
// BL depth for setting alpha0 outside BL
const dimensionedScalar BLdepth
(
    transferProps.lookupOrDefault<dimensionedScalar>
    (
        "BLdepth",
        dimensionedScalar("BLdepth", dimLength, scalar(0))
    )
);


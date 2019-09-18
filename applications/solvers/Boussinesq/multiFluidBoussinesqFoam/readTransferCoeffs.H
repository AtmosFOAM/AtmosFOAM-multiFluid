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
const dimensionedScalar cloudRadiusMax(transferProps.lookup("cloudRadiusMax"));
const dimensionedScalar cloudRadiusMin(transferProps.lookup("cloudRadiusMin"));
const scalar dragCoeff(readScalar(transferProps.lookup("dragCoeff")));

// Note if there are no transfers from any sources
const Switch noTransfers(transferProps.lookup("noTransfers"));

// Diffusion coefficient for inter-partition mass transfer
const dimensionedScalar Ksigma = noTransfers? dimlessZero : 
      dimensionedScalar(transferProps.lookup("Ksigma"));

// The amount of divergence transfer used (between zero and one)
const scalar divTransfer = noTransfers? scalar(0) :
      readScalar(transferProps.lookup("divTransfer"));
if (divTransfer < 0 || divTransfer > 1)
{
    FatalErrorIn("multiFluidBoussinesqFoam")
        << "divTransfer read in from file system/transferProperties should be "
        << "between zero and one, not " << divTransfer << exit(FatalError);
}

// Is mass transfer decided based on vertical velocity
const Switch wTransfer = noTransfers? Switch(false) : 
      transferProps.lookup("wTransfer");
const dimensionedScalar wTransferTimeScale = wTransfer ? 
      dimensionedScalar(transferProps.lookup("wTransferTimeScale")) :
      dimensionedScalar("transferProps", dimTime, scalar(0));

// Is mass transfer decided based on grad(b)
const Switch buoyancyTransfer = noTransfers? Switch(false) : 
      transferProps.lookup("buoyancyTransfer");
const dimensionedScalar bTransferDiffusivity
(
    transferProps.lookup("bTransferDiffusivity")
);

// Minimum sigma for each partition
const scalar minSigma = noTransfers? scalar(0) :
      readScalar(transferProps.lookup("minSigma"));

// Transfer to keep flagged cells uniform sigma
const Switch uniformSigma = noTransfers? Switch(false) :
      Switch(transferProps.lookup("uniformSigma"));

// Transfer buoyancy anomalies
const Switch transferAnomalies = noTransfers? Switch(false) :
        Switch(transferProps.lookup("transferAnomalies"));

// Bouyany of the fluid to transfer based on the laplacian
const Switch transferAnomaliesLaplacian = noTransfers? Switch(false) :
        Switch(transferProps.lookup("transferAnomaliesLaplacian"));
const dimensionedScalar bTransferTimeScale = transferAnomaliesLaplacian ?
      dimensionedScalar(transferProps.lookup("bTransferTimeScale")) :
      dimensionedScalar("", dimVelocity, scalar(0));
const dimensionedScalar bTransferMax = transferAnomaliesLaplacian ?
      dimensionedScalar(transferProps.lookup("bTransferMax")) :
      bZero;

const Switch transferAnomaliesb0 = noTransfers? Switch(false) :
        Switch(transferProps.lookup("transferAnomaliesb0"));
const Switch transferZeroAnomalies = noTransfers? Switch(false) :
        transferProps.lookupOrDefault<Switch>
        (
            "transferZeroAnomalies", Switch(false)
        );


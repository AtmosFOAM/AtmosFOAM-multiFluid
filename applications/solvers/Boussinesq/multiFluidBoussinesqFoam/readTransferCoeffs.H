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

// Note if there are no buoyancy anomalies transferred
const Switch noAnomalies(transferProps.lookup("noAnomalies"));

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
const dimensionedScalar wPrime = wTransfer ?
      dimensionedScalar(transferProps.lookup("wPrime")) :
      dimensionedScalar("dimensionedScalar", dimVelocity, scalar(0));
const dimensionedScalar btMax = wTransfer ?
      dimensionedScalar(transferProps.lookup("btMax")) : 
      dimensionedScalar("btMax", dimensionSet(0,1,-2,0,0), scalar(0));

// Is mass transfer decided based on grad(b)
const Switch buoyancyTransfer = noTransfers? Switch(false) :
      transferProps.lookup("buoyancyTransfer");
const dimensionedScalar bTransferDiffusivity =  buoyancyTransfer ?
      dimensionedScalar(transferProps.lookup("bTransferDiffusivity")) :
      dimensionedScalar("bTransferDiffusivity", dimless, scalar(0));

// Transfer w=0 (w^T = 0)
const Switch wTransfer0 = noTransfers? Switch(false) :
      transferProps.lookup("wTransfer0");

// Minimum sigma for each partition
const scalar minSigma = noTransfers? scalar(0) :
      readScalar(transferProps.lookup("minSigma"));

// Transfer to keep flagged cells uniform sigma
const Switch uniformSigma = noTransfers? Switch(false) :
      Switch(transferProps.lookup("uniformSigma"));

// Transfer buoyancy anomalies based on linear combination of mean fluid 
// buoyancies
const Switch transferAnomaliesLinComb = 
        noTransfers? Switch(false) : noAnomalies? Switch(false) : 
        transferProps.lookupOrDefault<Switch>
        (
            "transferAnomaliesLinComb", Switch(false)
        );
const dimensionedScalar bAnomalyLinearCoeff = transferAnomaliesLinComb ? 
      dimensionedScalar(transferProps.lookup("bAnomalyLinearCoeff")) :
      dimlessZero;
      
// Transfer buoyancy anomalies based on constant factor times mean buoyancy
const Switch transferAnomaliesMeanFactor = 
        noTransfers? Switch(false) : noAnomalies? Switch(false) : 
        transferProps.lookupOrDefault<Switch>
        (
            "transferAnomaliesMeanFactor", Switch(false)
        );
const scalar bAnomalyMeanCoeff = transferAnomaliesMeanFactor ? 
      readScalar(transferProps.lookup("bAnomalyMeanCoeff")) :
      scalar(0);

// Bouyany of the fluid to transfer based on the laplacian
const Switch transferAnomaliesLaplacian = 
        noTransfers? Switch(false) : noAnomalies? Switch(false) : 
        transferProps.lookupOrDefault<Switch>
        (
            "transferAnomaliesLaplacian", Switch(false)
        );
const dimensionedScalar bTransferTimeScale = transferAnomaliesLaplacian ?
      dimensionedScalar(transferProps.lookup("bTransferTimeScale")) :
      dimensionedScalar("", dimVelocity, scalar(0));
const dimensionedScalar bTransferMax = transferAnomaliesLaplacian ?
      dimensionedScalar(transferProps.lookup("bTransferMax")) :
      bZero;

const Switch transferAnomaliesb0 = 
        noTransfers? Switch(false) : noAnomalies? Switch(false) : 
        transferProps.lookupOrDefault<Switch>
        (
            "transferAnomaliesb0", Switch(false)
        );
const Switch transferZeroAnomalies = noTransfers? Switch(false) :
        transferProps.lookupOrDefault<Switch>
        (
            "transferZeroAnomalies", Switch(false)
        );

// Coefficient for entrainment rate proportional to 1/r
const scalar entrainmentCoefficient = transferProps.lookupOrDefault<scalar>
(
    "entrainmentCoefficient", scalar(0)
);


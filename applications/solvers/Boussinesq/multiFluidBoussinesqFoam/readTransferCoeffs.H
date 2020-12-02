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

// Transfer w=0 (w^T = 0) (scaled by sigma)
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
        noTransfers? Switch(false) :
        transferProps.lookupOrDefault<Switch>
        (
            "transferAnomaliesLinComb", Switch(false)
        );
const dimensionedScalar bAnomalyLinearCoeff = transferAnomaliesLinComb ? 
      dimensionedScalar(transferProps.lookup("bAnomalyLinearCoeff")) :
      dimlessZero;

// Transfer buoyancy anomalies based on constant factor times mean buoyancy
const Switch transferAnomaliesMeanFactor = 
        noTransfers? Switch(false) : 
        transferProps.lookupOrDefault<Switch>
        (
            "transferAnomaliesMeanFactor", Switch(false)
        );
const scalar bAnomalyMeanCoeff = transferAnomaliesMeanFactor ? 
      readScalar(transferProps.lookup("bAnomalyMeanCoeff")) :
      scalar(0);


Info << "\nReading environmentalProperties" << endl;

IOdictionary envProps
(
    IOobject
    (
        "environmentalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const dimensionedScalar nu(envProps.lookup("nu"));

const Switch adjustQ(envProps.lookupOrDefault<Switch>("adjustQ", false));

uniformDimensionedScalarField alpha
(
    IOobject("alpha", runTime.timeName(), mesh), 
    dimensionedScalar(envProps.lookup("alpha"))
);

uniformDimensionedScalarField gamma
(
    IOobject("gamma", runTime.timeName(), mesh), 
    dimensionedScalar(envProps.lookup("gamma"))
);

// Check that Ksigam is non-zero if alpha or nu are non-zero
if
(
    (mag(alpha.value() - Ksigma.value()) > SMALL)
 && (mag(Ksigma.value()) < SMALL)
)
{
    FatalErrorIn("multiFluidBoussinesqFoam if one out of alpha and Ksigma ")
         << "is zero then they must both be zero but alpha = "
         << alpha.value() << " and Ksigma = " << Ksigma.value()
         << exit(FatalError);
}
if
(
    (mag(nu.value() - Ksigma.value()) > SMALL)
 && (mag(Ksigma.value()) < SMALL)
)
{
    FatalErrorIn("multiFluidBoussinesqFoam if one out of nu and Ksigma ")
         << "is zero then they must both be zero but nu = "
         << nu.value() << " and Ksigma = " << Ksigma.value()
         << exit(FatalError);
}


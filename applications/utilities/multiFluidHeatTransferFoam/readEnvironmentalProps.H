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

// Viscosity
const dimensionedScalar nu(envProps.lookup("nu"));

// Thermal diffusivity
uniformDimensionedScalarField alpha
(
    IOobject("alpha", runTime.timeName(), mesh), 
    dimensionedScalar(envProps.lookup("alpha"))
);


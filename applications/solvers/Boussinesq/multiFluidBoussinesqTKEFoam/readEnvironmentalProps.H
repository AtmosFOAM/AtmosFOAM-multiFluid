Info << "\nReading environmentalProperties" << endl;

IOdictionary envProps
(
    IOobject
    (
        "environmentalProperties",
        runTime.system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Laminar viscosity
const dimensionedScalar nu(envProps.lookup("nu"));

// von-Karmen's constant
const dimensionedScalar vonKarmenKappa(envProps.lookup("vonKarmenKappa"));

// Maximum turbulent length scales
const List<dimensionedScalar> maxLength(envProps.lookup("maxLength"));


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

// Minimum and maximum turbulent length scales
const List<dimensionedScalar> minLength(envProps.lookup("minLength"));
const List<dimensionedScalar> maxLength(envProps.lookup("maxLength"));

// Constant to get TKE dissipation from k^1.5/l  (0.09 default)
const scalar Cmu = readScalar(envProps.lookup("Cmu"));

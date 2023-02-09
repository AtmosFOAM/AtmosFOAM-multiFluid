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

// Laminar viscosity and diffusivity
const dimensionedScalar nu(envProps.lookup("nu"));
const dimensionedScalar alpha(envProps.lookup("alpha"));

// von-Karmen's constant and roughness height
const dimensionedScalar vonKarmenKappa(envProps.lookup("vonKarmenKappa"));
const dimensionedScalar surfaceHeight(envProps.lookup("surfaceHeight"));

// Minimum and maximum turbulent length scales
const List<dimensionedScalar> minLength(envProps.lookup("minLength"));
const List<dimensionedScalar> maxLength(envProps.lookup("maxLength"));

// Constant to get TKE dissipation from k^1.5/l  (0.09 default)
const scalar Cmu = readScalar(envProps.lookup("Cmu"));

// Buoyancy TKE production coefficient (turbulent Prantl number)
const scalar Cb = readScalar(envProps.lookup("Cb"));

// Pressure gradient TKE production coefficient
const scalar Cp = readScalar(envProps.lookup("Cp"));


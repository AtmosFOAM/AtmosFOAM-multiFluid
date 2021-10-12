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
// Horizontal viscosity
const dimensionedScalar nuH
(
    envProps.lookupOrDefault<dimensionedScalar>("nuH", nu)
);

// Thermal diffusivity
uniformDimensionedScalarField alpha
(
    IOobject("alpha", runTime.timeName(), mesh), 
    dimensionedScalar(envProps.lookup("alpha"))
);
// Horizontal thermal diffusivity
const dimensionedScalar alphaH
(
    envProps.lookupOrDefault<dimensionedScalar>("alphaH", alpha)
);

// Surface Scalar fields for nu and alpha
surfaceScalarField nuf
(
    IOobject("nuf", runTime.timeName(), mesh), 
    magk*nu + (scalar(1)-magk)*nuH
);
surfaceScalarField alphaf
(
    IOobject("alphaf", runTime.timeName(), mesh), 
    magk*alpha + (scalar(1)-magk)*alphaH
);

// Scaling to go from distance to wall to thermal diffusivity of fluid 0
const dimensionedScalar alphaScale
(
    envProps.lookupOrDefault<dimensionedScalar>
    (
        "alphaScale", 
        dimensionedScalar("alphaScale", dimVelocity, scalar(0))
    )
);
// Maximum thermal diffusivity of fluid 0
const dimensionedScalar maxAlpha0
(
    envProps.lookupOrDefault<dimensionedScalar>
    (
        "maxAlpha0",
        dimensionedScalar("maxAlpha0", alpha.dimensions(), scalar(0))
    )
);


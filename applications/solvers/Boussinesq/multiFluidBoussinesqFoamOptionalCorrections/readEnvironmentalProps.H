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
const dimensionedScalar relaxbTime
(
    envProps.lookupOrDefault<dimensionedScalar>
    (
        "relaxBuoyancyTimeScale", 
        dimensionedScalar("relaxbTime", dimTime, scalar(0))
    )
);

uniformDimensionedScalarField alpha
(
    IOobject("alpha", runTime.timeName(), mesh), 
    dimensionedScalar(envProps.lookup("alpha"))
);

// Use diffusion of all fluids together or separate for each fluid
const Switch sharedDiffusion
(
    envProps.lookupOrDefault<Switch>("sharedDiffusion", false)
);

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

const dimensionedScalar Nsquared(envProps.lookup("Nsquared"));
const dimensionedScalar nu(envProps.lookup("nu"));
const dimensionedScalar alpha(envProps.lookup("alpha"));

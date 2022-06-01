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

// The partition names
const wordList partNames(transferProps.lookup("fluids"));
const label nParts = partNames.size();


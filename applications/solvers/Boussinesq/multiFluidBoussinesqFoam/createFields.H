// Read in prognostic variables
Info<< "Reading field P\n" << endl;
volScalarField P
(
    IOobject
    (
        "P",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// The partition names
const label nParts = 2;
wordList partNames(nParts);
partNames[0] = "stable";
partNames[1] = "buoyant";

Info << "Reading in sigma in each partition\n" << endl;
partitionedVolScalarField sigma
(
    "sigma", partNames, mesh, runTime.timeName()
);

// sigma on the faces
partitionedSurfaceScalarField sigmaf
(
    "sigmaf", partNames, fvc::interpolate(sigma[0])
);
for(label ip = 1; ip < nParts; ip++)
{
    sigmaf[ip] = fvc::interpolate(sigma[ip]);
}

Info << "Reading in b in each partition\n" << endl;
partitionedVolScalarField b
(
    "b", partNames, mesh, runTime.timeName(), sigma
);
b.mean().write();

Info << "Reading in cell centre velocity, u, in each partition\n" << endl;
partitionedVolVectorField u
(
    "u", partNames, mesh, runTime.timeName(), sigma
);
u.mean().write();

Info << "Calculating the flux in each partition\n"<<endl;
partitionedSurfaceScalarField volFlux
(
    "volFlux",
    partNames,
    surfaceScalarField
    (
        IOobject("volFlux"+'.'+partNames[0], runTime.timeName(), mesh),
        linearInterpolate(u[0]) & mesh.Sf()
    ),
    sigmaf
);
for(label ip = 1; ip < volFlux.size(); ip++)
{
    volFlux[ip] = linearInterpolate(u[ip]) & mesh.Sf();
}
volFlux.updateSum();
volFlux.mean().write();

Info << "Storing necessary old times and rates of change\n" << endl;
volFlux.storeTime();
b.storeTime();
sigma.storeTime();

// Surface scalar fields for hydrostatic boundary conditions
surfaceScalarField gradPcoeff
(
    IOobject("gradPcoeff", runTime.timeName(), mesh),
    fvc::interpolate(sigma.sum())
);
partitionedSurfaceScalarField bf
(
    "bf", partNames, fvc::interpolate(b[0], "b")*mesh.Sf().component(2), sigmaf
);
for(label ip = 1; ip < nParts; ip++)
{
    bf[ip] = fvc::interpolate(b[ip], "b")*mesh.Sf().component(2);
}
bf.updateSum();

Info << "Initialising the mass transfer between partitions\n" << endl;
volScalarField massTransferTmp
(
    IOobject("massTransfer", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    rateZero
);
PtrList<partitionedVolScalarField> massTransfer(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    massTransfer.set
    (
        jp, 
        new partitionedVolScalarField
        (
            "massTransfer."+partNames[jp], partNames, massTransferTmp
        )
    );
}

// Transfer terms for other variables (assuming nParts == 2)
partitionedVolScalarField transfer
(
    "transfer",
    partNames,
    volScalarField
    (
        IOobject("transfer", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("", dimless, scalar(0))
    )
);

partitionedSurfaceScalarField transferf
(
    "transferf",
    partNames,
    surfaceScalarField
    (
        IOobject("transfer", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("", dimless, scalar(0))
    )
);

// transferLocation = 1 where sigma is to be kept uniform
partitionedVolScalarField transferLocation
(
    "transferLocation",
    partNames,
    volScalarField
    (
        IOobject("transferLocation", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("", dimless, scalar(0))
    )
);

if (uniformSigma)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        transferLocation[ip] == volScalarField
        (
            IOobject("transferLocation."+partNames[ip], runTime.timeName(), mesh,
                    IOobject::MUST_READ, IOobject::NO_WRITE),
            mesh
        );
    }
}

Info << "Initialising the drag between partitions" << endl;
volScalarField cloudRadius
(
    "cloudRadius",
    max(cloudRadiusMin, cloudRadiusMax*sigma[0]*sigma[1])
);
volScalarField dragCommon = -dragCoeff/cloudRadius* mag(u[0] - u[1]);

// Rate of change of b (eg for radiative cooling)
uniformDimensionedScalarField Q
(
    IOobject("Q", runTime.timeName(), mesh, 
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    dimensionedScalar("Q", dimensionSet(0,1,-3,0,0), scalar(0))
);
Q.write();

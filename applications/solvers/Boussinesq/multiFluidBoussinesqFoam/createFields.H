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
sigma.sum().write();

// sigma on the faces
partitionedSurfaceScalarField sigmaf
(
    "sigmaf", partNames, linearInterpolate(sigma[0])
);
for(label ip = 1; ip < nParts; ip++)
{
    sigmaf[ip] = linearInterpolate(sigma[ip]);
}
sigmaf.updateSum();

// Pressure in each partition
partitionedVolScalarField Pi("Pi", partNames, P, sigma);

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

// fixed flux boundary conditions
wordList fluxBCs(mesh.boundaryMesh().size(), "calculated");
forAll(mesh.boundaryMesh(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    if (btype == "wall" || btype == "symmetryPlane" || btype == "empty")
    {
        fluxBCs[patchi] = "fixedValue";
    }
}

Info << "Calculating the flux in each partition\n"<<endl;
partitionedSurfaceScalarField volFlux
(
    "volFlux",
    partNames,
    surfaceScalarField
    (
        IOobject("volFlux"+'.'+partNames[0], runTime.timeName(), mesh),
        linearInterpolate(u[0]) & mesh.Sf(),
        fluxBCs
    ),
    sigmaf
);
for(label ip = 1; ip < volFlux.size(); ip++)
{
    volFlux[ip] = linearInterpolate(u[ip]) & mesh.Sf();
}
volFlux.updateSum();
volFlux.sum().write();

Info << "Storing necessary old times and rates of change\n" << endl;
volFlux.storeTime();
b.storeTime();
sigma.storeTime();

// Divergence (for post processing)
partitionedVolScalarField divu
(
    "divu", partNames, fvc::div(volFlux[0]), sigma
);
for(label ip = 1; ip < nParts; ip++)
{
    divu[ip] = fvc::div(volFlux[ip]);
}
divu.write();

// Velocity on the faces
partitionedSurfaceVectorField Uf
(
    "Uf", partNames, linearInterpolate(u[0]), sigmaf
);
for(label ip = 0; ip < nParts; ip++)
{
    Uf[ip] = linearInterpolate(fvc::reconstruct(volFlux[ip]));
    Uf[ip] += (volFlux[ip] - (Uf[ip] & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());
}
Uf.updateSum();
Uf.write();

// Surface scalar fields for hydrostatic boundary conditions
const surfaceScalarField gradPcoeff
(
    IOobject("gradPcoeff", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("1", dimless, scalar(1))
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
const volScalarField transferLocation
(
    IOobject("transferLocation", runTime.constant(), mesh,
             IOobject::READ_IF_PRESENT, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("", dimless, scalar(0))
);

if (uniformSigma)
{
    // Check that transferLocation is non-zero somewhere
    if (abs(max(transferLocation).value()) < SMALL)
    {
        FatalErrorIn("multiFluidBoussinesqFoam") << " transferProperty uniformSigma set to true but no transferLocation set to non-zero"
            << exit(FatalError);
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
             IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    dimensionedScalar("Q", dimensionSet(0,1,-3,0,0), scalar(0))
);
Q.write();

// Buoyancy difference of the transferred fluid
partitionedVolScalarField bt
(
    "bt",
    partNames,
    volScalarField
    (
        IOobject("bt", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("", b[0].dimensions(), scalar(0)),
        "fixedValue"
    ),
    sigma
);

// Volume flux difference of the transferred fluid
partitionedSurfaceScalarField volFluxt
(
    "volFluxt",
    partNames,
    surfaceScalarField
    (
        IOobject("volFluxt", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("", volFlux[0].dimensions(), scalar(0)),
        fluxBCs
    ),
    sigmaf
);


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

Info << "Reading in sigma in each partition\n" << endl;
partitionedVolScalarField sigma
(
    "sigma", partNames, mesh, runTime.timeName()
);
sigma.sum().write();

Info << "Setting sigma on the faces\n" << endl;
partitionedSurfaceScalarField sigmaf
(
    "sigmaf", partNames, linearInterpolate(sigma[0]), IOobject::NO_WRITE
);
interpolate(sigmaf, sigma);

Info << "Reading in b in each partition\n" << endl;
partitionedVolScalarField b
(
    "b", partNames, mesh, runTime.timeName(), sigma
);
b.mean().write();

// Surface scalar fields for hydrostatic boundary conditions
const surfaceScalarField gradPcoeff
(
    IOobject("gradPcoeff", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("1", dimless, scalar(1))
);
partitionedSurfaceScalarField bf
(
    "bf", partNames, fvc::interpolate(b[0], "b")*mesh.Sf().component(2),
    sigmaf, IOobject::NO_WRITE
);
for(label ip = 1; ip < nParts; ip++)
{
    bf[ip] = fvc::interpolate(b[ip], "b")*mesh.Sf().component(2);
}
bf.updateSum();

Info << "Reading pressure anomaly in each partition\n" << endl;
partitionedVolScalarField p
(
    "p",
    partNames,
    mesh,
    runTime.timeName(),
    sigma
);

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

Info << "Reading/calculating the volFlux in each partition\n"<<endl;
partitionedSurfaceScalarField volFlux
(
    "volFlux",
    partNames,
    surfaceScalarField
    (
        IOobject("volFlux"+'.'+partNames[0], runTime.timeName(), mesh,
                IOobject::READ_IF_PRESENT),
        linearInterpolate(u[0]) & mesh.Sf(),
        fluxBCs //wordList(mesh.boundaryMesh().size(), "calculated")
    ),
    sigmaf
);
for(label ip = 1; ip < nParts; ip++)
{
    IOobject volFluxIOobject
    (
        volFlux[ip].name(), runTime.timeName(), mesh, IOobject::MUST_READ
    );
    if (volFluxIOobject.headerOk())
    {
        volFlux[ip] = surfaceScalarField(volFluxIOobject, mesh);
    }
    else
    {
        volFlux[ip] = linearInterpolate(u[ip]) & mesh.Sf();
    }
}
volFlux.updateSum();

Info << "Storing necessary old times and rates of change\n" << endl;
volFlux.storeTime();
u.storeTime();
b.storeTime();
sigma.storeTime();

volScalarTransferField massTransfer
(
    "massTransfer",
    partNames,
    volScalarField
    (
        IOobject("massTransfer", runTime.timeName(), mesh),
        mesh,
        rateZero,
        "fixedValue"
    )
);
massTransfer.write();

// Mass transfers on faces
surfaceScalarTransferField massTransferf
(
    "massTransferf",
    partNames,
    surfaceScalarField
    (
        IOobject("massTransferf", runTime.timeName(), mesh),
        mesh,
        rateZero
    ),
    IOobject::NO_WRITE
);

// Drag between fluids
surfaceScalarTransferField drag
(
    "drag",
    partNames,
    surfaceScalarField
    (
        IOobject("drag", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("", dimensionSet(0,3,-2,0,0), scalar(0))
    ),
    IOobject::NO_WRITE
);

// Rate of change of b (eg for radiative cooling)
volScalarField Q
(
    IOobject("Q", runTime.timeName(), mesh, 
             IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("Q", dimensionSet(0,1,-3,0,0), scalar(0))
);
// Total initial heat in the domain
dimensionedScalar totalHeat0
(
    "totalHeat0", dimensionSet(0,1,-2,0,0), scalar(0)
);
for(label ip = 0; ip < nParts; ip++)
{
    totalHeat0 += fvc::domainIntegrate(sigma[ip]*b[ip])/sum(mesh.V());
}
Info << "totalHeat0 = " << totalHeat0 << endl;

// Buoyancy of the transferred fluid
volScalarTransferField bT
(
    "bT",
    partNames,
    volScalarField
    (
        IOobject("bT", runTime.timeName(), mesh, IOobject::NO_READ),
        mesh,
        dimensionedScalar("", b[0].dimensions(), scalar(0))
    )/*,
    IOobject::NO_WRITE*/
);

// Velocity anomaly of the transferred fluid
surfaceScalarTransferField volFluxT
(
    "volFluxT",
    partNames,
    surfaceScalarField
    (
        IOobject("volFluxT", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("", volFlux[0].dimensions(), scalar(0)),
        fluxBCs
    ),
    IOobject::NO_WRITE
);

// Cell centred velocity of the transferred fluid (for output)
//volVectorTransferField uT("uT", partNames, u[0]);

// Compact Laplacian of the velocity
surfaceScalarField laplacianU
    = (fvc::interpolate(fvc::laplacian(nuf,u.updateSum())) & mesh.Sf())
    + mesh.magSf()*(nuf*fvc::snGrad
    (
        fvc::div
        (
            volFlux.updateSum()
         - (linearInterpolate(u.sum()) & mesh.Sf())
        )
    ));

// Divergence of each fluid
partitionedVolScalarField divu
(
    "divu",
    partNames,
    fvc::div(volFlux.mean()),
    sigma,
    IOobject::NO_WRITE
);

// Distance to a wall for each cell
wallDist wd(mesh);
const volScalarField wallDistance("wallDistance", max(wd.y() - BLdepth, lengthZero));

// Additional eddy diffusivity of heat for fluid 0
const surfaceScalarField alpha0
(
    "alpha0",
    min(linearInterpolate(wallDistance)*alphaScale, maxAlpha0)
);


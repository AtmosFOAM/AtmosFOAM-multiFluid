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
        fluxBCs
    ),
    sigmaf
);
for(label ip = 1; ip < nParts; ip++)
{
    IOobject volFluxIOobject
    (
        volFlux[ip].name(), runTime.timeName(), mesh, IOobject::MUST_READ
    );
    if (volFluxIOobject.typeHeaderOk<surfaceScalarField>(false))
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
    ),
    IOobject::NO_WRITE
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
volVectorTransferField uT("uT", partNames, u[0]);

// Divergence of each fluid
partitionedVolScalarField divu
(
    "divu",
    partNames,
    fvc::div(volFlux.mean()),
    sigma,
    IOobject::NO_WRITE
);

// Turbulent kinetic energy (TKE) in each fluid
partitionedVolScalarField TKE
(
    "TKE", partNames, mesh, runTime.timeName(), sigma
);
TKE.mean().write();
TKE.storeTime();

// TKE of the transferred fluid
volScalarTransferField tkeT
(
    "tkeT", partNames,
    volScalarField
    (
        IOobject("tkeT", runTime.timeName(), mesh, IOobject::NO_READ),
        mesh,
        dimensionedScalar("", TKE[0].dimensions(), scalar(0))
    ),
    IOobject::NO_WRITE
);
volScalarTransferField tkeAdd
(
    "tkeAdd", partNames,
    volScalarField
    (
        IOobject("tkeAdd", runTime.timeName(), mesh, IOobject::NO_READ),
        mesh,
        dimensionedScalar("", TKE[0].dimensions(), scalar(0))
    )
);
volScalarTransferField tkeRemove
(
    "tkeRemove", partNames,
    volScalarField
    (
        IOobject("tkeRemove", runTime.timeName(), mesh, IOobject::NO_READ),
        mesh,
        dimensionedScalar("", TKE[0].dimensions(), scalar(0))
    )
);

// Fixed turbulent length scale based on distance to a wall
wallDist wd(mesh);
partitionedVolScalarField turbLength
(
    "turbLength", partNames, wd.y()*vonKarmenKappa, IOobject::NO_WRITE
);
for(label ip = 0; ip < nParts; ip++)
{
    turbLength[ip] = min(turbLength[ip], maxLength[ip]);
}
turbLength.write();

// Turbulent viscosity in each fluid
partitionedSurfaceScalarField nuTurb
(
    "nuTurb",
    partNames,
    surfaceScalarField
    (
        IOobject("nuTurb", runTime.timeName(), mesh),
        mesh,
        nu
    ),
    sigmaf
);
for(label ip = 0; ip < nParts; ip++)
{
    nuTurb[ip] = nu + linearInterpolate(turbLength[ip]*sqrt(TKE[ip]));
}
nuTurb.write();

// Compact Laplacian of the velocity
surfaceScalarField laplacianU
    = (fvc::interpolate(fvc::laplacian(nuTurb[0],u.updateSum())) & mesh.Sf())
    + mesh.magSf()*nuTurb[0]*fvc::snGrad
    (
        fvc::div
        (
            volFlux.updateSum()
         - (linearInterpolate(u.sum()) & mesh.Sf())
        )
    );


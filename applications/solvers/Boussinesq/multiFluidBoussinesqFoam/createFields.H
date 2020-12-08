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

Info << "Setting sigma on the faces\n" << endl;
partitionedSurfaceScalarField sigmaf
(
    "sigmaf", partNames, linearInterpolate(sigma[0])//, IOobject::NO_WRITE
);
for(label ip = 1; ip < nParts; ip++)
{
    sigmaf[ip] = linearInterpolate(sigma[ip]);
}
sigmaf.updateSum();

Info << "Setting pressure anomaly in each partition\n" << endl;
partitionedVolScalarField p
(
    "p",
    partNames,
    volScalarField
    (
        IOobject("p", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("", P.dimensions(), scalar(0)),
        "zeroGradient"
    ),
    sigma
);

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
volFlux.mean().write();

Info << "Storing necessary old times and rates of change\n" << endl;
volFlux.storeTime();
u.storeTime();
b.storeTime();
sigma.storeTime();

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
    "bf", partNames, fvc::interpolate(b[0], "b")*mesh.Sf().component(2),
    sigmaf, IOobject::NO_WRITE
);
for(label ip = 1; ip < nParts; ip++)
{
    bf[ip] = fvc::interpolate(b[ip], "b")*mesh.Sf().component(2);
}
bf.updateSum();

// The mass transfer specific to diffusion
partitionedVolScalarField massTransferD
(
    "massTransferD",
    partNames, 
    volScalarField
    (
        IOobject("massTransferD", runTime.timeName(), mesh),
        mesh,
        rateZero,
        "fixedValue"
    )//,    IOobject::NO_WRITE
);

PtrList<partitionedVolScalarField> massTransfer(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    massTransfer.set
    (
        jp, 
        new partitionedVolScalarField
        (
            "massTransfer."+partNames[jp], partNames, massTransferD[0]
        )
    );
    massTransfer[jp].write();
}

// Transfer terms for other variables (assuming nParts == 2)
// transfer[0] = dt*massTransfer[0][1]/sigma[1]
partitionedVolScalarField transfer
(
    "transfer",
    partNames,
    volScalarField
    (
        IOobject("transfer", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("", dimless, scalar(0))
    ),
    IOobject::NO_WRITE
);

partitionedVolScalarField transferD
(
    "transferD",
    partNames,
    volScalarField
    (
        IOobject("transferD", runTime.timeName(), mesh),
        mesh,
        rateZero
    ),
    IOobject::NO_WRITE
);

// Transfer terms for other variables (assuming nParts == 2)
//S[0] = dt*massTransfer[0][1]/sigma[0]
partitionedVolScalarField S
(
    "S",
    partNames,
    volScalarField
    (
        IOobject("S", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("", dimless, scalar(0))
    ),
    IOobject::NO_WRITE
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

// Write out initial massTransfer for uniformSigma
if (uniformSigma)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        sigma.ddt()[ip] = -fvc::div(volFlux[ip], sigma[ip], "sigma");
    }

    massTransfer[0][1] = transferLocation*max(sigma.ddt()[0], rateZero)
                       + (1-transferLocation)*massTransfer[0][1];
    massTransfer[1][0] = transferLocation*max(sigma.ddt()[1], rateZero)
                       + (1-transferLocation)*massTransfer[1][0];
    massTransfer[0].write();
    massTransfer[1].write();
}


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

// Buoyancy difference of the transferred fluid
partitionedVolScalarField bt
(
    "bt",
    partNames,
    volScalarField
    (
        IOobject("bt", runTime.timeName(), mesh, 
                IOobject::NO_READ, IOobject::AUTO_WRITE),
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
    sigmaf,
    IOobject::NO_WRITE
);

// Transfer terms on the faces
partitionedSurfaceScalarField massTransferf
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
partitionedSurfaceScalarField transferf
(
    "transferf",
    partNames,
    surfaceScalarField
    (
        IOobject("transferf", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("", dimless, scalar(0))
    ),
    IOobject::NO_WRITE
);
partitionedSurfaceScalarField Sf
(
    "Sf",
    partNames,
    surfaceScalarField
    (
        IOobject("Sf", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("", dimless, scalar(0))
    ),
    IOobject::NO_WRITE
);
partitionedSurfaceScalarField transferDf
(
    "transferDf",
    partNames,
    surfaceScalarField
    (
        IOobject("transferDf", runTime.timeName(), mesh),
        mesh,
        rateZero
    ),
    IOobject::NO_WRITE
);

partitionedVolScalarField divu
(
    "divu", partNames, fvc::div(sigmaf[0]*volFlux[0])
);
divu[1] = fvc::div(sigmaf[1]*volFlux[1]);
divu.write();

// Compact Laplacian of the velocity
surfaceScalarField laplacianU
    = (fvc::interpolate(fvc::laplacian(u.updateSum())) & mesh.Sf())
    + mesh.magSf()*(fvc::snGrad
    (
        fvc::div(volFlux.updateSum()
      - (linearInterpolate(u.sum()) & mesh.Sf()))
    ));


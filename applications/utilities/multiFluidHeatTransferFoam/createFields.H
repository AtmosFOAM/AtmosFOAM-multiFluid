// Read in fields

Info << "Reading in sigma in each partition\n" << endl;
partitionedVolScalarField sigma
(
    "sigma", partNames, mesh, runTime.timeName()
);

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

Info << "Reading/calculating the volFlux in each partition\n"<<endl;
partitionedSurfaceScalarField volFlux
(
    "volFlux",
    partNames,
    surfaceScalarField
    (
        IOobject("volFlux"+'.'+partNames[0], runTime.timeName(), mesh,
                IOobject::READ_IF_PRESENT),
        mesh,
        dimensionedScalar("", dimensionSet(0,3,-1,0,0), scalar(0)),
        wordList(mesh.boundaryMesh().size(), "calculated")
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
}
volFlux.updateSum();

// Read in turbulent viscosity if it exists
partitionedSurfaceScalarField nuTurb
(
    "nuTurb",
    partNames,
    surfaceScalarField
    (
        IOobject("nuTurb"+'.'+partNames[0], runTime.timeName(), mesh,
                IOobject::READ_IF_PRESENT),
        mesh,
        alpha,
        wordList(mesh.boundaryMesh().size(), "calculated")
    ),
    sigmaf
);
for(label ip = 1; ip < nParts; ip++)
{
    IOobject nuTurbIOobject
    (
        nuTurb[ip].name(), runTime.timeName(), mesh, IOobject::MUST_READ
    );
    if (nuTurbIOobject.typeHeaderOk<surfaceScalarField>(false))
    {
        nuTurb[ip] = surfaceScalarField(nuTurbIOobject, mesh);
    }
}

// Create heat transfer in each fluid
partitionedSurfaceScalarField heatTransferf
(
    "heatTransferf", partNames, -alpha*fvc::snGrad(b[0]), sigmaf,
    IOobject::NO_WRITE
);

partitionedVolVectorField heatTransfer
(
    "heatTransfer", partNames, fvc::reconstruct(heatTransferf[0]*mesh.magSf()),
    sigma, IOobject::NO_WRITE
);




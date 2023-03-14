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
        IOobject("volFlux."+partNames[0], runTime.timeName(), mesh,
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
    if (volFluxIOobject.headerOk())
    {
        volFlux[ip] = surfaceScalarField(volFluxIOobject, mesh);
    }
}
volFlux.updateSum();

Info << "Reading in turbulent diffusivity if it exists\n" << endl;
partitionedSurfaceScalarField alphaTurb
(
    "alphaTurb",
    partNames,
    surfaceScalarField
    (
        IOobject("alphaTurb."+partNames[0], runTime.timeName(), mesh,
                IOobject::READ_IF_PRESENT),
        mesh,
        alpha,
        wordList(mesh.boundaryMesh().size(), "calculated")
    ),
    sigmaf
);
for(label ip = 1; ip < nParts; ip++)
{
    IOobject alphaTurbIOobject
    (
        alphaTurb[ip].name(), runTime.timeName(), mesh, IOobject::MUST_READ
    );
    if (alphaTurbIOobject.headerOk())
    {
        alphaTurb[ip] = surfaceScalarField(alphaTurbIOobject, mesh);
    }
}

Info << "Creating heat transfer in each fluid\n" << endl;
partitionedSurfaceScalarField heatTransferf
(
    "heatTransferf", partNames, -alpha*fvc::snGrad(b[0])*sigmaf[0]
);

partitionedVolVectorField heatTransfer
(
    "heatTransfer", partNames, fvc::reconstruct(heatTransferf[0]*mesh.magSf())
);




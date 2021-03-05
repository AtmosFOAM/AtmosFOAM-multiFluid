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
    "sigmaf", partNames, linearInterpolate(sigma[0]), IOobject::NO_WRITE
);
for(label ip = 1; ip < nParts; ip++)
{
    sigmaf[ip] = linearInterpolate(sigma[ip]);
}
sigmaf.updateSum();

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

/*// Laplacian of pressure based on random numbers in order to aid convergence
Random rand(1);
volScalarField h("h", fvc::laplacian(p[0]));
forAll(h, cellI)
{
    h[cellI] = rand.scalarAB(-1,1);
}
h -= fvc::domainIntegrate(h)/dimensionedScalar("", dimVol,gSum(mesh.V()));
h *= dt*fvc::domainIntegrate(fvc::div(bf.sum()))
    /dimensionedScalar("", dimVol*h.dimensions(), scalar(1));

// Dummy pressure to balance p
volScalarField ph("ph", p[0]);
{
    fvScalarMatrix pEqn(fvm::laplacian(ph) - h/dt);
    pEqn.setReference(0,0);
    pEqn.solve();
}
ph.write();
h.write();
*/

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

PtrList<partitionedVolScalarField> massTransfer(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    massTransfer.set
    (
        jp, 
        new partitionedVolScalarField
        (
            "massTransfer."+partNames[jp],
            partNames,
            volScalarField
            (
                IOobject("massTransfer", runTime.timeName(), mesh),
                mesh,
                rateZero,
                "fixedValue"
            )
        )
    );
    massTransfer[jp].write();
}

// Diffusive transfers
PtrList<partitionedVolScalarField> massTransferD(nParts);
for(label jp = 0; jp < nParts; jp++)
{
    massTransferD.set
    (
        jp, 
        new partitionedVolScalarField
        (
            "massTransferD."+partNames[jp],
            partNames,
            volScalarField
            (
                IOobject("massTransferD", runTime.timeName(), mesh),
                mesh,
                rateZero,
                "fixedValue"
            ),
            IOobject::NO_WRITE
        )
    );
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

// Buoyancy of the transferred fluid
partitionedVolScalarField bT
(
    "bT",
    partNames,
    volScalarField
    (
        IOobject("bT", runTime.timeName(), mesh, 
                IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("", b[0].dimensions(), scalar(0)),
        "fixedValue"
    ),
    sigma
);

// Velocity anomaly of the transferred fluid
partitionedSurfaceScalarField volFluxT
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
    sigmaf,
    IOobject::NO_WRITE
);

// Compact Laplacian of the velocity
surfaceScalarField laplacianU
    = (fvc::interpolate(fvc::laplacian(u.updateSum())) & mesh.Sf())
    + mesh.magSf()*(fvc::snGrad
    (
        fvc::div(volFlux.updateSum()
      - (linearInterpolate(u.sum()) & mesh.Sf()))
    ));


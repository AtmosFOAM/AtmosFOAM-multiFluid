volScalarField ExnerPrime
(
    IOobject("ExnerPrime", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("zero", dimless, scalar(0)),
    "zeroGradient"
);

// The partition names
const label nParts = 2;
wordList partNames(nParts);
partNames[0] = "stable";
partNames[1] = "buoyant";

Info<< "Reading field Exner in each partition\n" << endl;
partitionedVolScalarField Exner
(
    "Exner", partNames, mesh, runTime.timeName()
);

Info << "Reading in sigma in each partition\n" << endl;
partitionedVolScalarField sigma
(
    "sigma", partNames, mesh, runTime.timeName()
);

Info << "Setting rho for each partition\n" << endl;
partitionedVolScalarField rho
(
    "rho",
    partNames,
    volScalarField
    (
        IOobject("rho", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("rho", dimDensity, scalar(1)),
        sigma[0].boundaryField().types()
    ),
    sigma
);
partitionedVolScalarField sigmaRho = rho.timesSigma();

Info << "Reading in theta in each partition\n" << endl;
partitionedVolScalarField theta
(
    "theta", partNames, mesh, runTime.timeName(), sigmaRho
);

Info << "Calculating density and density in each partition and density weighted theta sum\n" << endl;
for(label ip = 0; ip < partNames.size(); ip++)
{
    rho[ip] = pRef/(R*theta[ip])*pow(Exner[ip], (1-kappa)/kappa);
}
rho.updateSum();
sigmaRho = rho.timesSigma();
sigmaRho.updateSum();
theta.updateSum();
sigmaRho.write();

Info << "Initialising sigmaRhof and sigmaf - values on the faces\n" << endl;
partitionedSurfaceScalarField sigmaRhof
(
    "sigmaRhof", partNames, fvc::interpolate(sigmaRho[0], "sigmaRho")
);
for(label ip = 1; ip < partNames.size(); ip++)
{
    sigmaRhof[ip] = fvc::interpolate(sigmaRho[ip], "sigmaRho");
}
sigmaRhof.updateSum();

partitionedSurfaceScalarField thetaf
(
    "thetaf", partNames, fvc::interpolate(theta[0], "theta"), sigmaRhof
);
for(label ip = 1; ip < partNames.size(); ip++)
{
    thetaf[ip] = fvc::interpolate(theta[ip], "theta");
}
thetaf.updateSum();

partitionedSurfaceScalarField sigmaf
(
    "sigmaf", partNames, fvc::interpolate(sigma[0], "sigmaRho")
);
for(label ip = 1; ip < partNames.size(); ip++)
{
    sigmaf[ip] = fvc::interpolate(sigma[ip], "sigmaRho");
}
sigmaf.updateSum();

Info << "Reading in the velocity in each partition\n" << endl;
partitionedSurfaceVectorField Uf
(
    "Uf", partNames, mesh, runTime.timeName(), sigmaRhof
);

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

Info << "Calculating the mass flux and volume flux in each partition\n"<<endl;
partitionedSurfaceScalarField volFlux
(
    "volFlux",
    partNames,
    surfaceScalarField
    (
        IOobject("volFlux"+'.'+partNames[0], runTime.timeName(), mesh),
        Uf[0] & mesh.Sf(),
        fluxBCs
    ),
    sigmaRhof
);
for(label ip = 1; ip < volFlux.size(); ip++)
{
    volFlux[ip] = Uf[ip] & mesh.Sf();
}
volFlux.updateSum();

Info << "Reading in BCs and calculating the cell centre velocity, u, in each partition\n" << endl;
partitionedVolVectorField u
(
    "u", partNames, mesh, runTime.timeName(), sigmaRho
);
for(label ip = 0; ip < partNames.size(); ip++)
{
    u[ip] = fvc::reconstruct(volFlux[ip]);
}

Info << "Storing necessary old times and rates of change\n" << endl;
volFlux.storeTime();
theta.storeTime();
sigmaRho.storeTime();
rho.storeTime();
sigma.storeTime();

Info << "Initialising the drag between partitions" << endl;
surfaceScalarField cloudRadius
(
    "cloudRadius",
    max(cloudRadiusMin, cloudRadiusMax*sigmaf[0]*sigmaf[1])
);
surfaceScalarField dragCommon = -dragCoeff*sigmaRhof.sum()/cloudRadius
                         * mag(Uf[0] - Uf[1]);


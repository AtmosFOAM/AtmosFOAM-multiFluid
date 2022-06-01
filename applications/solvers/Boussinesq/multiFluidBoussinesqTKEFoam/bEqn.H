// flux form b equation without diffusion or transfers

// Diffusion common for all fluids
volScalarField bDiffusion = sigma[0]*fvc::laplacian(nuTurb[0], b[0]);
for(label ip = 1; ip < nParts; ip++)
{
    bDiffusion += sigma[ip]*fvc::laplacian(nuTurb[ip], b[ip]);
}

for(label ip = 0; ip < nParts; ip++)
{
    b.ddt()[ip] = Q - fvc::div(volFlux[ip], b[ip], "b")
                + b[ip]*fvc::div(volFlux[ip])
                + bDiffusion;

    // Update b according to rate of change
    b[ip] = b[ip].oldTime() + dt*
    (
        (1-offCentre)*b.ddt()[ip].oldTime()
      + offCentre*b.ddt()[ip]
    );

    b[ip].correctBoundaryConditions();

    // bouyancy force on faces
    bf[ip] = fvc::interpolate(b[ip], "b")*mesh.Sf().component(2);
}

b.updateSum();
bf.updateSum();



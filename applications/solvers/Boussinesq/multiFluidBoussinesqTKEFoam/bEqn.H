// flux form b equation without transfers

for(label ip = 0; ip < nParts; ip++)
{
    b.ddt()[ip] = Q - fvc::div(volFlux[ip], b[ip], "b")
                + b[ip]*fvc::div(volFlux[ip])
                + Cb*fvc::laplacian(nuTurb[ip], b[ip]);

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


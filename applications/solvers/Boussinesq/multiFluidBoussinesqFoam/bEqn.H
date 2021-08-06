// flux form b equation without diffusion or transfers
for(label ip = 0; ip < nParts; ip++)
{
    b.ddt()[ip] = Q - fvc::div(volFlux[ip], b[ip], "b")
                + b[ip]*fvc::div(volFlux[ip])
                + fvc::laplacian(alphaf,b[ip]);

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

//Info << "b[1] - b[0] goes from " << min(b[1] - b[0]).value() << " to "
//     << max(b[1] - b[0]).value() << endl;



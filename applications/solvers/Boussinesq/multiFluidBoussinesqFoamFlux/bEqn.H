// flux form b equation without diffusion or transfers
for(label ip = 0; ip < nParts; ip++)
{
    b.ddt()[ip] = Q - fvc::div(flux[ip], b[ip], "b");

    // Update b according to rate of change and transfers
    b[ip] = sigma[ip].oldTime()*b[ip].oldTime() + dt*
    (
        (1-offCentre)*b.ddt()[ip].oldTime()
      + offCentre*b.ddt()[ip]
    );

    // Transfers
    for(label jp=0; jp < nParts; jp++) if (ip != jp)
    {
        b[ip] += dt*(massTransfer[jp][ip]*bT[jp] - massTransfer[ip][jp]*bT[ip]);
    }

    // Update b
    b[ip] /= max(sigma[ip], minSigmaDiv);
    b[ip].correctBoundaryConditions();

    // bouyancy force on faces
    bf[ip] = fvc::interpolate(b[ip], "b")*mesh.Sf().component(2);
}

b.updateSum();
bf.updateSum();

Info << "b[1] goes from " << min(b[1]).value() <<  " to "
    << max(b[1]).value() << endl;

volScalarField bDiff("bDiff", b[1] - b[0]);
bDiff.write();
Info << "bDiff goes from " << min(bDiff).value() << " to "
     << max(bDiff).value() << endl;


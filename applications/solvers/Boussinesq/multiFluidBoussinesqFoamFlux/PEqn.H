// Pressure equation
for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix PEqn
    (
        fvc::div(flux.sum())
      + offCentre*dt*fvc::div(bf.sum())
      - offCentre*dt*fvm::laplacian(P)
    );
    PEqn.setReference(0,0);
    PEqn.solve();
}

// Back-substitutions
for(label ip = 0; ip < nParts; ip++)
{
    flux.ddt()[ip] += sigmaf[ip]*(bf[ip] - fvc::snGrad(P)*mesh.magSf());
    flux[ip] += offCentre*dt*sigmaf[ip]*(bf[ip] - fvc::snGrad(P)*mesh.magSf());
}
flux.updateSum();

partitionedVolVectorField u2
(
    "u2",
    partNames,
    fvc::reconstruct(flux[0]/sigmaf[0]),
    sigma
);
u2[1] = fvc::reconstruct(flux[1]/sigmaf[1]);
u2.write();
volVectorField u2Diff("u2Diff2", u2[1] - u2[0]);
u2Diff.write();

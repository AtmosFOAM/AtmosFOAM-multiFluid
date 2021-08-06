// Pressure equation
for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix PEqn
    (
        fvc::div(volFlux.sum())
      + offCentre*dt*fvc::div(bf.sum())
      - offCentre*dt*fvm::laplacian(P)
    );
    PEqn.setReference(mesh.nCells()/2, scalar(0));
    PEqn.solve();
}

// Back-substitutions
for(label ip = 0; ip < nParts; ip++)
{
    volFlux.ddt()[ip] += bf[ip] - fvc::snGrad(P)*mesh.magSf();
    volFlux[ip] += offCentre*dt*(bf[ip] - fvc::snGrad(P)*mesh.magSf());
}

/*volScalarField divFlux("divFlux", fvc::div(volFlux.sum()));
volScalarField divb("divb", offCentre*dt*fvc::div(bf.sum()));
volScalarField del2P("del2P", offCentre*dt*fvc::laplacian(P));
divFlux.write();
divb.write();
del2P.write();

volFlux.updateSum();
divFlux.rename("divFlux2");
divFlux = fvc::div(volFlux.sum());
divFlux.write();
*/

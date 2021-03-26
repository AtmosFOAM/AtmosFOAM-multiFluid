// Pressure equation in each fluid
for(label ip = 0; ip < nParts; ip++)
{
    //sumM[ip] = rateZero;
    //sumM.transferMass(massTransfer, scalar(1));
    sumM[ip] = massTransfer(1-ip,ip) - massTransfer(ip,1-ip);
}
sumM.updateSum();

for(label ip = 0; ip < nParts; ip++)
{
    volFlux.ddt()[ip] -= bf[ip] - fvc::snGrad(P[ip])*mesh.magSf();
    volFlux[ip] -= offCentre*dt*
    (
        bf[ip] - fvc::snGrad(P[ip])*mesh.magSf()
    );

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix PEqn
        (
           // fvm::Sp(1/gamma, P[ip])
            fvc::ddt(sigma[ip])
          + fvc::div(volFlux[ip])
          + offCentre*dt*fvc::div(bf[ip])
          - offCentre*dt*fvm::laplacian(P[ip])
         == sumM[ip]/max(sigma[ip], dimlessSMALL)
        );
        PEqn.setReference(label(mesh.nCells()/2),scalar(0));
        PEqn.solve();
    }

    // Back-substitutions
    volFlux.ddt()[ip] += bf[ip] - fvc::snGrad(P[ip])*mesh.magSf();
    volFlux[ip] += offCentre*dt*(bf[ip] - fvc::snGrad(P[ip])*mesh.magSf());
}
volFlux.updateSum();
// Update velocities based on the volFlux
for(label ip = 0; ip < nParts; ip++)
{
    u[ip] = fvc::reconstruct(volFlux[ip]);
}
u.updateSum();
P.updateSum();

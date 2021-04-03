// Pressure anomalies in each fluid
// Pi = P + pi
// where pi = -gamma div(u_i) + gamma sum_j(sigma_j div(u_j)

/*for(label ip = 0; ip < nParts; ip++)
{
    //sumM[ip] = rateZero;
    //sumM.transferMass(massTransfer, scalar(1));
    sumM[ip] = massTransfer(1-ip,ip) - massTransfer(ip,1-ip);
}
sumM.updateSum();
*/
if (gamma.value() > SMALL)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        volFlux[ip] += offCentre*dt*fvc::snGrad(p[ip])*mesh.magSf();
        volFlux.ddt()[ip] += fvc::snGrad(p[ip])*mesh.magSf();

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                fvm::Sp(1/gamma, p[ip])
              + fvc::div(volFlux[ip])
              - offCentre*dt*fvm::laplacian(p[ip])
              //- offCentre*dt*fvc::laplacian(p.sum())
              //== sumM[ip]/max(sigma[ip], dimlessSMALL)
            );
            pEqn.solve();
        }
    }
    p.updateSum();
    Info << "p.sum goes from " << min(p.sum()).value() << " to "
         << max(p.sum()).value() << endl;
    for(label ip = 0; ip < nParts; ip++)
    {
        //p[ip] -= p.sum();
        // Back-substitutions
        volFlux[ip] -= offCentre*dt*fvc::snGrad(p[ip])*mesh.magSf();
        volFlux.ddt()[ip] -= fvc::snGrad(p[ip])*mesh.magSf();
    }
}

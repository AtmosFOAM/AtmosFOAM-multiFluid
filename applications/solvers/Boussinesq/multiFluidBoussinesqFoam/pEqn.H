// Pressure anomalies in each fluid
// Pi = P + pi
// where pi = -gamma div(u_i) + gamma sum_j(sigma_j div(u_j)

if (gamma.value() > SMALL)
{
    // Iterate to achieve p.sum = 0
    for(label ipIt = 0; ipIt < 1; ipIt++)
    {
        for(label ip = 0; ip < nParts; ip++)
        {
            volFlux[ip] += offCentre*dt*(fvc::snGrad(p[ip])*mesh.magSf());
            divu[ip] = fvc::div(volFlux[ip]);
        }
        divu.updateSum();

        for(label ip = 0; ip < nParts; ip++)
        {
            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::Sp(1/gamma, p[ip])
                  - divu.sum()
                  + fvc::div(volFlux[ip])
                  - offCentre*dt*fvm::laplacian(p[ip])
                  //+ offCentre*dt*fvc::laplacian(p[ip])
                );
                pEqn.solve();
            }
            volFlux[ip] -= offCentre*dt*(fvc::snGrad(p[ip])*mesh.magSf());
        }
        p.updateSum();
        Info << "p.sum goes from " << min(p.sum()).value() << " to "
             << max(p.sum()).value() << endl;
/*        for(label ip = 0; ip < nParts; ip++)
        {
            p[ip] -= p.sum();
        }
*/    }
}


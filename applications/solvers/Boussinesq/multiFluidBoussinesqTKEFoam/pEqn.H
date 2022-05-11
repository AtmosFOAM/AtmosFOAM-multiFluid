// Pressure anomalies in each fluid
// Pi = P + pi
// where pi = -gamma div(u_i) + gamma sum_j(sigma_j div(u_j)

if (gamma.value() > SMALL)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        divu[ip] = fvc::div(volFlux[ip]);
        volFlux[ip] += offCentre*dt*fvc::snGrad(p[ip])*mesh.magSf();
        volFlux.ddt()[ip] += fvc::snGrad(p[ip])*mesh.magSf();
    }
    divu.updateSum();
    
    for(label ip = 0; ip < nParts; ip++)
    {
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                fvm::Sp(1/gamma, p[ip])
              + fvc::div(volFlux[ip])
              - offCentre*dt*fvm::laplacian(p[ip])
              - divu.sum()
            );
            pEqn.solve();
        }

        // Back-substitutions
        volFlux[ip] -= offCentre*dt*fvc::snGrad(p[ip])*mesh.magSf();
        volFlux.ddt()[ip] -= fvc::snGrad(p[ip])*mesh.magSf();
    }
    
    p.updateSum();
    Info << "p.sum goes from " << min(p.sum()).value() << " to "
         << max(p.sum()).value() << endl;
}

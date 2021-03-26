// Pressure anomalies in each fluid
// Pi = P + pi
// where pi = -gamma div(u_i) + gamma sum_j(sigma_j div(u_j)

if (gamma.value() > SMALL)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        volFlux[ip] += offCentre*dt*fvc::snGrad(p[ip])*mesh.magSf();
    }
    
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
            );
            pEqn.solve();
        }
        volFlux[ip] -= offCentre*dt*fvc::snGrad(p[ip])*mesh.magSf();
        p.updateSum();
        Info << "p.sum goes from " << min(p.sum()).value() << " to "
             << max(p.sum()).value() << endl;
    }
}


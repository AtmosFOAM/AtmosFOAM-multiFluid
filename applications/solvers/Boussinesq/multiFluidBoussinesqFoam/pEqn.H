// Pressure anomalies in each fluid
// Pi = P + pi
// where pi = -gamma div(u_i) + gamma sum_j(sigma_j div(u_j)

if (gamma.value() > VSMALL)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                fvm::Sp(1/gamma, p[ip])
              + fvc::div(volFlux[ip])
              + offCentre*dt*fvc::laplacian(p[ip])
              - offCentre*dt*fvm::laplacian(p[ip])
            );

            pEqn.solve();
        }
    }
    p.updateSum();
    for(label ip = 0; ip < nParts; ip++)
    {
        p[ip] -= p.sum();
    }
}


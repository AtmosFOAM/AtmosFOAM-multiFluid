// Pressure anomalies in each fluid
// Pi = P + pi
// where pi = -gamma div(u_i) + gamma sum_j(sigma_j div(u_j)

if (gamma.value() > SMALL)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        flux.ddt()[ip] += sigmaf[ip]*fvc::snGrad(p[ip])*mesh.magSf();
    }

    // Iterate as each p is dependent on the divergence in both fluids
    for(int pCorr = 0; pCorr < 1; pCorr++)
    {
        partitionedVolScalarField divu
        (
            "divu", partNames, fvc::div(flux[0]/sigmaf[0]), sigma
        );
        divu[1] = fvc::div(flux[1]/sigmaf[1]);
        divu.updateSum();
        divu.write();
    
        for(label ip = 0; ip < nParts; ip++)
        {
            flux[ip] += sigmaf[ip]*offCentre*dt*fvc::snGrad(p[ip])*mesh.magSf();
            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::Sp(1/gamma, p[ip])
                  + fvc::div(flux[ip]/sigmaf[ip])
                  - offCentre*dt*fvm::laplacian(p[ip])
                  - divu.sum()
                );

                pEqn.solve();
            }
            flux[ip] -= sigmaf[ip]*offCentre*dt*fvc::snGrad(p[ip])*mesh.magSf();
        }
        p.updateSum();
        Info << "p.sum goes from " << min(p.sum()).value() << " to "
             << max(p.sum()).value() << endl;
    }
    for(label ip = 0; ip < nParts; ip++)
    {
        flux.ddt()[ip] -= sigmaf[ip]*fvc::snGrad(p[ip])*mesh.magSf();
    }
}


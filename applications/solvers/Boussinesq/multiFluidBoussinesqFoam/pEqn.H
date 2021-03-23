// Pressure anomalies in each fluid
// Pi = P + pi
// where pi = -gamma div(u_i) + gamma sum_j(sigma_j div(u_j)

if (mag(gamma.value()) > SMALL)
{
    const dimensionedScalar oneByGamma = gamma.value() > SMALL? 1/gamma :
            dimensionedScalar("", dimensionSet(0,-2,1,0,0), scalar(0));

    // Iterate to achieve p.sum = 0
    for(label ipIt = 0; ipIt < 1; ipIt++)
    {
        for(label ip = 0; ip < nParts; ip++)
        {
            volFlux[ip] -= offCentre*dt*
            (
                bf[ip] - fvc::snGrad(p[ip])*mesh.magSf()
            );
            divu[ip] = fvc::div(volFlux[ip]);
        }
        divu.updateSum();

        for(label ip = 0; ip < nParts; ip++)
        {
            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::Sp(oneByGamma, p[ip])
                  //- divu.sum()
                  + fvc::div(volFlux[ip])
                  - offCentre*dt*fvm::laplacian(p[ip])
                );
                if (oneByGamma.value() < SMALL)
                {
                    pEqn.setReference(int(mesh.nCells()/2),scalar(0));
                }
                pEqn.solve();
            }
            volFlux[ip] += offCentre*dt*
            (
                bf[ip] - fvc::snGrad(p[ip])*mesh.magSf()
            );
        }
        p.updateSum();
        Info << "p.sum goes from " << min(p.sum()).value() << " to "
             << max(p.sum()).value() << endl;
    }
}


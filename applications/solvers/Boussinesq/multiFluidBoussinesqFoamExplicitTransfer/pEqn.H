// Pressure anomalies in each fluid
// Pi = P + pi
// where pi = -gamma div(u_i) + gamma sum_j(sigma_j div(u_j)

/*Info << "volFlux[1] - volFlux[0] goes from "
     << min(volFlux[1] - volFlux[0]).value() << " to "
     << max(volFlux[1] - volFlux[0]).value() << endl;
*/
if (gamma.value() > SMALL)
{
    // Iterate to achieve p.sum = 0
    for(label ip = 0; ip < nParts; ip++)
    {
        volFlux[ip] -= offCentre*dt*(fvc::snGrad(p[ip])*mesh.magSf());
    }
    for(label ipIt = 0; ipIt < 1; ipIt++)
    {
        volScalarField divuSum = sigma[0]*fvc::div(volFlux[0])
                               + sigma[1]*fvc::div(volFlux[1]);

        for(label ip = 0; ip < nParts; ip++)
        {
            volFlux[ip] += offCentre*dt*(fvc::snGrad(p[ip])*mesh.magSf());

//            p[ip] == p[ip] + ph;
            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::Sp(1/gamma, p[ip])
//                   + offCentre*h - ph/gamma
                  - divuSum
                  + fvc::div(volFlux[ip])
                  - offCentre*dt*fvm::laplacian(p[ip])
                );
                pEqn.solve();
            }
//            p[ip] == p[ip] - ph;
            volFlux[ip] -= offCentre*dt*(fvc::snGrad(p[ip])*mesh.magSf());
        }
        p.updateSum();
        Info << "p.sum goes from " << min(p.sum()).value() << " to "
             << max(p.sum()).value() << endl;
    }
    for(label ip = 0; ip < nParts; ip++)
    {
        volFlux.ddt()[ip] -= fvc::snGrad(p[ip])*mesh.magSf();
    }
}


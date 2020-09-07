// Separate pressure equation for each fluid

if (gamma.value() > VSMALL)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            // Correction to previous pressure as
            // Pi = -gamma div(u_i) + gamma sum(sigma_i div(u_i)
            fvScalarMatrix PEqn
            (
                fvm::Sp(1/gamma, Pi[ip])
              - Pi.sum()/gamma
              + fvc::div(volFlux[ip])
              + offCentre*dt*fvc::laplacian(Pi[ip])
              - offCentre*dt*fvm::laplacian(Pi[ip])
            );
            
            PEqn.solve();
        }
    }
    Pi.updateSum();
}


// Separate pressure equation for each fluid

if (gamma.value() > VSMALL)
{
    // Correction to the previous fluid pressure as
    // sigma Pi = - gamma div(sigma u_i)
    for(label ip = 0; ip < nParts; ip++)
    {
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
/*            fvScalarMatrix PEqn
            (
                fvm::Sp(sigma[ip]/gamma, Pi[ip])
              + fvc::div(sigmaf[ip]*volFlux[ip])
              + offCentre*dt*fvc::laplacian(sigmaf[ip],Pi[ip])
              - offCentre*dt*fvm::laplacian(sigmaf[ip],Pi[ip])
            );*/
            fvScalarMatrix PEqn
            (
                fvm::Sp(1/gamma, Pi[ip])
              + fvc::div(volFlux[ip])
              + offCentre*dt*fvc::laplacian(Pi[ip])
              - offCentre*dt*fvm::laplacian(Pi[ip])
            );
            PEqn.solve();
        }
    }
    Pi.updateSum();
}


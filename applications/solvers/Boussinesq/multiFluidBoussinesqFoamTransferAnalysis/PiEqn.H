// Separate pressure equation for each fluid

if (gamma.value() > VSMALL)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            // Correction to the previous fluid pressure as
            //  Pi = - gamma div(sigma u_i)
            // This version conserves the energy E where E_i = K_i + Pi
            // and E = sum(sigma_i E_i)
            /*fvScalarMatrix PEqn
            (
                fvm::Sp(1/gamma, Pi[ip])
              + fvc::div(sigmaf[ip]*(volFlux[ip]-volFlux.mean()))
              //+ fvc::div(sigmaf[ip]*volFlux[ip])
              + offCentre*dt*fvc::laplacian(sigmaf[ip],Pi[ip])
              - offCentre*dt*fvm::laplacian(sigmaf[ip],Pi[ip])
            );*/
            // Correction to previous pressure as
            // Pi = -gamma div(u_i)
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


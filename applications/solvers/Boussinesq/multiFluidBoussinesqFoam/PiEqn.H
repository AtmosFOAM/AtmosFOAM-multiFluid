// Separate pressure equation for each fluid

if (gamma.value() > VSMALL)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix PEqn
            (
                fvm::Sp(1/gamma,Pi[ip])
              + fvc::div(volFlux[ip])
              + offCentre*dt*fvc::div(bf[ip])
              - offCentre*dt*fvm::laplacian(Pi[ip])
    //            fvm::ddt(Pi[ip])
    //          + gamma*fvc::div(volFlux[ip])
    //          - gamma*offCentre*dt*fvm::laplacian(Pi[ip])
            );
            PEqn.setReference(0,0);
            PEqn.solve();
        }
        
        // Back substitutions
        volFlux[ip] -= offCentre*dt*fvc::snGrad(Pi[ip])*mesh.magSf();
    }
}

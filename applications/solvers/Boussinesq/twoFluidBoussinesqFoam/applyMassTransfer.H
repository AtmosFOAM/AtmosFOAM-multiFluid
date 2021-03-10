// Apply mass transfer terms (operator split) to sigma equation
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            sigma[ip] += dt*
            (
                massTransfer[jp][ip] - massTransfer[ip][jp]
            );
        }
    }
    sigma[ip].correctBoundaryConditions();
}
sigma.updateSum();

// Apply mass transfer terms (operator split) to sigmaf
for(label ip = 0; ip < nParts; ip++)
{
    sigmaf[ip] += dt*(massTransferf[1-ip] - massTransferf[ip]);
}
sigmaf.updateSum();
volFlux.updateSum();
bf.updateSum();
b.updateSum();


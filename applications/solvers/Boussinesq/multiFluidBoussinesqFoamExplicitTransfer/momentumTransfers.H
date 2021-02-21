// Velocity of the fluid transferred
volFluxT[0] = volFlux[0];
volFluxT[1] = volFlux[1];

if (!noTransfers && wTransfer0)
{
    // Velocity anomalies of fluid transferred
    volFluxT[0] += sigmaf[0]*wPrimeCoeff*mag(volFlux[0]*mesh.Sf().component(2))
                   *mesh.Sf().component(2)/sqr(mesh.magSf());
    volFluxT[1] -= sigmaf[1]*wPrimeCoeff*mag(volFlux[1]*mesh.Sf().component(2))
                   *mesh.Sf().component(2)/sqr(mesh.magSf());
}

// Transfer the velocity
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp = 0; jp < nParts; jp++) if (ip != jp)
    {
        volFlux[ip].oldTime() += dt/max(sigmaf[ip].oldTime(), minSigmaDiv)*
        (
            fvc::interpolate(massTransfer[jp][ip])
                *(volFluxT[jp] - volFlux[ip].oldTime())
          - fvc::interpolate(massTransfer[ip][jp])
            *(volFluxT[ip] - volFlux[ip].oldTime())
        );
    }
}


volFluxT[0] = flux[0]/max(sigmaf[0], minSigmaDiv);
volFluxT[1] = flux[1]/max(sigmaf[1], minSigmaDiv);

if (!noTransfers && wTransfer0)
{
    // Velocity difference of fluid transferred
    volFluxT[0] += wPrimeCoeff*mag(flux[0]*mesh.Sf().component(2))
                   *mesh.Sf().component(2)/sqr(mesh.magSf());
    volFluxT[1] -= wPrimeCoeff*mag(flux[1]*mesh.Sf().component(2))
                   *mesh.Sf().component(2)/sqr(mesh.magSf());
}


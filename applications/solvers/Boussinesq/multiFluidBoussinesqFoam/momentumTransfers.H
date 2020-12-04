if (!noTransfers)
{
    // Transfers due to diffusion coefficients for momentum and sigma being 
    // different. Apply these transfers independently
    if (mag(nu - Ksigma).value() > SMALL)
    {
        surfaceScalarField A = dt*(nu - Ksigma)/Ksigma
                            *fvc::interpolate(massTransferD[0], "transfer");

        volFlux[0] += (volFlux[0] - A*volFlux[1])/(1 - A);

        A = dt*(nu - Ksigma)/Ksigma
            *fvc::interpolate(massTransferD[1], "transfer");

        volFlux[1] += (volFlux[1] - A*volFlux[0])/(1 - A);
    }

    // Velocity difference of fluid transferred
    if (wTransfer0) for(label ip = 0; ip < nParts; ip++)
    {
        volFluxt[ip] = -sigmaf[ip]*volFlux[ip]
                       *mesh.Sf().component(2)/mesh.magSf();
    }

    // Implicit updates of volFlux for each fluid separately
    volFlux[0] =
    (
        (1 + transferf[0])*
        (
            volFlux[0] - Sf[0]*volFluxt[0] + transferf[1]*volFluxt[1]
        )
      + transferf[1]*
        (
            volFlux[1] + transferf[0]*volFluxt[0] - Sf[1]*volFluxt[1]
        )
    )/
    (1 + transferf.sum());

    volFlux[1] = 
    (
        volFlux[1]
      + (transferf[0])*volFlux[0]
      + transferf[0]*volFluxt[0]
      - Sf[1]*volFluxt[1]
    )
    /(1 + transferf[0]);
}


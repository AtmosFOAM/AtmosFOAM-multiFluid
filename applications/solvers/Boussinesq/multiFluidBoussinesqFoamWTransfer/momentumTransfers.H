if (!noTransfers || dragCoeff > 0)
{
    // Velocity difference of fluid transferred
    for(label ip = 0; ip < nParts; ip++)
    {
        volFluxt[ip] = -volFlux[ip]*mesh.Sf().component(2)/mesh.magSf();
    }

    // Implicit updates of volFlux for each fluid separately
    volFlux[0] =
    (
        (1 + transferf[0] + dt*dragCommon*sigmaf[0])*
        (
            volFlux[0] - Sf[0]*volFluxt[0] + transferf[1]*volFluxt[1]
        )
      + (transferf[1] + dt*dragCommon*sigmaf[1])*
        (
            volFlux[1] + transferf[0]*volFluxt[0] - Sf[1]*volFluxt[1]
        )
    )/
    (1 + transferf.sum() + dt*dragCommon);

    volFlux[1] = 
    (
        volFlux[1]
      + (transferf[0] + dt*sigmaf[0]*dragCommon)*volFlux[0]
      + transferf[0]*volFluxt[0]
      - Sf[1]*volFluxt[1]
    )
    /(1 + transferf[0] + dt*sigmaf[0]*dragCommon);
}


if (!noTransfers || dragCoeff > 0)
{
    // Velocity difference of fluid transferred
    if (wTransfer0) for(label ip = 0; ip < nParts; ip++)
    {
//        volFluxt[ip] = -sigmaf[ip]*wPrimeCoeff*volFlux[ip]
//                       *mesh.Sf().component(2)/mesh.magSf();
        volFluxt[0] =  sigmaf[0]*wPrimeCoeff*
                        mag(volFlux[0]*mesh.Sf().component(2))
                       *mesh.Sf().component(2)/sqr(mesh.magSf());
        volFluxt[1] =  -sigmaf[1]*wPrimeCoeff*
                        mag(volFlux[1]*mesh.Sf().component(2))
                       *mesh.Sf().component(2)/sqr(mesh.magSf());
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


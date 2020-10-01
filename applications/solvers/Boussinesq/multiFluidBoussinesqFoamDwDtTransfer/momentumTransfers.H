if (!noTransfers || dragCoeff > 0)
{
    partitionedSurfaceScalarField volFluxPreTransfer = volFlux;

    // Velocity difference of fluid transferred
    for(label ip = 0; ip < nParts; ip++)
    {
        if (wTransferType == "transferMean")
        {
            volFluxt[ip] *= 0;
        }
        else if (wTransferType == "transferMeanCell")
        {
            volFluxt[ip] = (volFlux.mean()-volFlux[ip])*mesh.Sf().component(2)/mesh.magSf();
        }
        else if (wTransferType == "transferZero")
        {
            volFluxt[ip] = -volFlux[ip]*mesh.Sf().component(2)/mesh.magSf();
        }
        else if (wTransferType == "transferPdf")
        {
            volFluxt[ip] = (volFluxT[ip]-volFlux[ip])*mesh.Sf().component(2)/mesh.magSf();
        }
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

    
    // Implicit scheme for transferring the mean of a cell.
    // surfaceScalarField transferDiff = dt*(massTransferf[0]-massTransferf[1]);

    // volFlux[0] = 
    // (
    //     sigmaf[0]*sigmaf[1]*volFluxPreTransfer[0]
    //   + (sqr(sigmaf[0])*volFluxPreTransfer[0] - sqr(sigmaf[1])*volFluxPreTransfer[1])*transferDiff
    //   - (sigmaf[0]*volFluxPreTransfer[0] + sigmaf[1]*volFluxPreTransfer[1])*sqr(transferDiff)
    // ) / 
    // (
    //     sigmaf[0]*sigmaf[1]
    //   + (sqr(sigmaf[0]) - sqr(sigmaf[1]))*transferDiff
    //   - sqr(transferDiff)
    // );

    // volFlux[1] = 
    // (
    //     sigmaf[0]*sigmaf[1]*volFluxPreTransfer[1]
    //   + (sqr(sigmaf[0])*volFluxPreTransfer[0] - sqr(sigmaf[1])*volFluxPreTransfer[1])*transferDiff
    //   - (sigmaf[0]*volFluxPreTransfer[0] + sigmaf[1]*volFluxPreTransfer[1])*sqr(transferDiff)
    // ) / 
    // (
    //     sigmaf[0]*sigmaf[1]
    //   + (sqr(sigmaf[0]) - sqr(sigmaf[1]))*transferDiff
    //   - sqr(transferDiff)
    // );
}


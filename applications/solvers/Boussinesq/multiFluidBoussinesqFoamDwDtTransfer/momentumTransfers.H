if (!noTransfers || dragCoeff > 0)
{
    partitionedSurfaceScalarField volFluxPreTransfer = volFlux;

    // Velocity difference of fluid transferred
    for(label ip = 0; ip < nParts; ip++)
    {
        volFluxt[ip] = -volFlux[ip]*mesh.Sf().component(2)/mesh.magSf();
        // volFluxt[ip] = -sigmaf[ip]*volFlux[ip]*mesh.Sf().component(2)/mesh.magSf();
        // volFluxt[ip] = (volFlux.mean()-volFlux[ip])*mesh.Sf().component(2)/mesh.magSf();
        // volFluxt[ip] = sigmaf[ip]*sigmaf[1-ip]*(volFlux[1-ip]-volFlux[ip])*mesh.Sf().component(2)/mesh.magSf();
        //volFluxt[ip] = (volFluxT[ip]-volFlux[ip])*mesh.Sf().component(2)/mesh.magSf();
        // volFluxt[ip] *= 0;
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

    // volFlux[1] = 
    // (
    //     volFlux[1]
    //   + (transferf[0] + dt*sigmaf[0]*dragCommon)*volFlux[0]
    //   + transferf[0]*volFluxt[0]
    //   - Sf[1]*volFluxt[1]
    // )
    // /(1 + transferf[0] + dt*sigmaf[0]*dragCommon);

    volFlux[1] =
    (
        (1 + transferf[1] + dt*dragCommon*sigmaf[1])*
        (
            volFlux[1] - Sf[1]*volFluxt[1] + transferf[0]*volFluxt[0]
        )
      + (transferf[0] + dt*dragCommon*sigmaf[0])*
        (
            volFlux[0] + transferf[1]*volFluxt[1] - Sf[0]*volFluxt[0]
        )
    )/
    (1 + transferf.sum() + dt*dragCommon);

    
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

    forAll (volFlux[0], celli)
    {
        // if (1)
        if (mag(volFlux[0][celli]-volFluxPreTransfer[0][celli]) > 1e-8)
        {
            Info << endl;
            Info << "stage transfers " << endl;
            Info << "t=" << runTime.timeName() << ", ucorr=" << ucorr << ", corr=" << corr << " " << endl;
            Info << celli << endl;
            Info << "sigma0: " << sigmaf[0][celli] << endl;
            Info << "sigma1: " << sigmaf[1][celli] << endl;
            Info << "S0: " << Sf[0][celli] << endl;
            Info << "S1: " << Sf[1][celli] << endl;
            Info << "w0: " << volFluxPreTransfer[0][celli]/2e7 << " to " << volFlux[0][celli]/2e7 << endl;
            Info << "w1: " << volFluxPreTransfer[1][celli]/2e7 << " to " << volFlux[1][celli]/2e7 << endl;
            Info << "w.ddt0: " << volFlux.ddt()[0][celli]/2e7 << endl;
            Info << "w.ddt1: " << volFlux.ddt()[1][celli]/2e7 << endl;
        }
    }
}


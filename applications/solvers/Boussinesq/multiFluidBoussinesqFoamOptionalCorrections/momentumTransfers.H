if (!noTransfers || dragCoeff > 0)
{
    // Velocity difference of fluid transferred
    for(label ip = 0; ip < nParts; ip++)
    {
        volFluxt[ip] = -volFlux[ip]*mesh.Sf().component(2)/mesh.magSf();
    }
    
/*    if (runTime.timeName() == "136")
    {
        forAll(volFlux[0], celli)
        {
            Info << "Test1: " << (1 + transferf[0][celli] + 2*dragCommon[celli]*sigmaf[0][celli]) << endl;
            Info << "Test2: " << (volFlux[0][celli] - Sf[0][celli]*volFluxt[0][celli] + transferf[1][celli]*volFluxt[1][celli]) << endl;
            
            Info << "Test12: " << ((1 + transferf[0][celli] + 2*dragCommon[celli]*sigmaf[0][celli])*(volFlux[0][celli] - Sf[0][celli]*volFluxt[0][celli] + transferf[1][celli]*volFluxt[1][celli])) << endl;
            
            Info << "Test3: " << (transferf[1][celli] + 2*dragCommon[celli]*sigmaf[1][celli]) << endl;
            Info << "Test4: " << (volFlux[1][celli] + transferf[0][celli]*volFluxt[0][celli] - Sf[1][celli]*volFluxt[1][celli]) << endl;
            
            Info << "Test3a: " << (transferf[1][celli]) << endl;
            Info << "Test3b: " << (2*dragCommon[celli]) << endl;
            Info << "Test3c: " << (sigmaf[1][celli]) << endl;
            
            Info << "Test34: " << ((transferf[1][celli] + 2*dragCommon[celli]*sigmaf[1][celli])*(volFlux[1][celli] + transferf[0][celli]*volFluxt[0][celli] - Sf[1][celli]*volFluxt[1][celli])) << endl;
            
            Info << "Test1234: " << ((1 + transferf[0][celli] + 2*dragCommon[celli]*sigmaf[0][celli])*(volFlux[0][celli] - Sf[0][celli]*volFluxt[0][celli] + transferf[1][celli]*volFluxt[1][celli])+(transferf[1][celli] + 2*dragCommon[celli]*sigmaf[1][celli])*(volFlux[1][celli] + transferf[0][celli]*volFluxt[0][celli] - Sf[1][celli]*volFluxt[1][celli])) << endl;
            
            Info << "Test5: " << (1 + transferf.sum()[celli] + 2*dragCommon[celli]) << endl;
            Info << "FinalTest: " << (((1 + transferf[0][celli] + 2*dragCommon[celli]*sigmaf[0][celli])*(volFlux[0][celli] - Sf[0][celli]*volFluxt[0][celli] + transferf[1][celli]*volFluxt[1][celli]) + (transferf[1][celli] + 2*dragCommon[celli]*sigmaf[1][celli])*(volFlux[1][celli] + transferf[0][celli]*volFluxt[0][celli] - Sf[1][celli]*volFluxt[1][celli]))/(1 + transferf.sum()[celli] + 2*dragCommon[celli])) << endl;
        }
    }
    Info << "here1" << endl;*/
    
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


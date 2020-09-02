// The transfer terms are only formulated for nParts == 2
if (nParts != 2)
{
    FatalErrorIn("massTransfer.H")
        << " mass transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}

// Set all mass transfers to zero
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=ip+1; jp < nParts; jp++)
    {
        massTransfer[ip][jp] == rateZero;
        massTransfer[jp][ip] == rateZero;
    }
}


/*
massTransferTmp = 0.5*Ksigma*fvc::laplacian
(
    sigma[1] -sigma[0]
);
massTransferredLaplacianSigma[0] = dt*max(massTransferTmp, rateZero);
massTransferredLaplacianSigma[1] = dt*min(massTransferTmp, rateZero);*/


massTransferredCloudRadius[0] = dt*max
(
    sigma[0]*u[1].component(2)*entrainmentCoefficient/cloudRadiusMax,
    rateZero
);
massTransferredCloudRadius[1] = dt*max
(
    sigma[1]*u[0].component(2)*entrainmentCoefficient/cloudRadiusMax,
    rateZero
);


massTransferredDivTransfer[0] = dt*sigma[0]*max(-fvc::div(volFlux[0]), rateZero);
massTransferredDivTransfer[1] = dt*sigma[1]*max(-fvc::div(volFlux[1]), rateZero);



/*
massTransferredLaplacianBuoyancy[0] = dt*sigma[0]*bTransferDiffusivity*min
(
    fvc::laplacian(b[0])/max
    (
        b[0],
        dimensionedScalar("", dimensionSet(0,1,-2,0,0), SMALL)
    ),
    dimensionedScalar("", dimensionSet(0,-2,0,0,0), scalar(0))
);
massTransferredLaplacianBuoyancy[1] = dt*sigma[1]*bTransferDiffusivity*max
(
    fvc::laplacian(b[1])/max
    (
        b[1],
        dimensionedScalar("", dimensionSet(0,1,-2,0,0), SMALL)
    ),
    dimensionedScalar("", dimensionSet(0,-2,0,0,0), scalar(0))
);*/


Info << 1 << endl;
massTransferredBuoyancyTransfer[0] = dt*max
(
    1e-6*cloudRadiusMax/dt * b[0]/max
    (
        magSqr(u[0]),
        dimensionedScalar("", dimensionSet(0,2,-2,0,0), SMALL)
    ),
    rateZero
);
massTransferredBuoyancyTransfer[1] = dt*max
(
    1e-6*cloudRadiusMax/dt * b[1]/max
    (
        magSqr(u[1]),
        dimensionedScalar("", dimensionSet(0,2,-2,0,0), SMALL)
    ),
    rateZero
);


/*Info << 2 << endl;
for(label ip = 0; ip < nParts; ip++)
{
    bMinMax[ip] = -(fvc::reconstruct
    (
        volFlux.ddt()[ip]
      + volFlux[ip]*(Sf[ip]/dt + 1/wTransferTimeScale)
    ) & vector(0,0,1));
}
massTransferredWTransfer[0] = dt*max
(
    0*massTransfer[0][1],
    sigma[0]/max(bt[0], bSMALL)*
    (
        -bMinMax[0]/wTransferTimeScale + b.ddt()[0]
    )
);
massTransferredWTransfer[1] = dt*max
(
    0*massTransfer[1][0],
    sigma[1]/min(bt[1], -bSMALL)*
    (
        -bMinMax[1]/wTransferTimeScale + b.ddt()[1]
    )
);*/



volVectorField DuDt0 = fvc::reconstruct(volFlux.ddt()[0] + (mesh.Sf() &
(
   fvc::interpolate
   (
       fvc::div(volFlux[0], u[0], "u")
     - u[0]*fvc::div(volFlux[0]),
       "convection"
   )
)) );
volVectorField DuDt1 = fvc::reconstruct(volFlux.ddt()[1] + (mesh.Sf() &
(
   fvc::interpolate
   (
       fvc::div(volFlux[1], u[1], "u")
     - u[1]*fvc::div(volFlux[1]),
       "convection"
   )
)));

massTransferredDwDt[0] = dt*DuDt0.component(2)/dimensionedScalar("", dimensionSet(0,1,-1,0,0), scalar(1));
massTransferredDwDt[1] = dt*DuDt1.component(2)/dimensionedScalar("", dimensionSet(0,1,-1,0,0), scalar(1));


massTransferredDwDtTransfer[0] += dt*sigma[0] * max
(
   - DuDt0.component(2)/
     min
     (
         u[0].component(2),
         dimensionedScalar("", dimensionSet(0,1,-1,0,0), -SMALL)
     ),
     rateZero
);
massTransferredDwDtTransfer[1] += dt*sigma[1] * max
(
   - DuDt1.component(2)/
     max
     (
         u[0].component(2),
         dimensionedScalar("", dimensionSet(0,1,-1,0,0), SMALL)
     ),
     rateZero
);



massTransferredDwDtPdfTransfer1[0] = sigma[0] * 0.5 *
(
    1 + Foam::erf
    (
        (u[0].component(2) + dt*DuDt0.component(2)) /
        max
        (
            1.41 * mag(u[0].component(2)-u.mean().component(2)),
            dimensionedScalar("", dimensionSet(0,1,-1,0,0), SMALL)
        )
    )
);
massTransferredDwDtPdfTransfer1[1] = sigma[1] * 0.5 *
(
    1 - Foam::erf
    (
        (u[1].component(2) - dt*DuDt1.component(2)) /
        max
        (
            1.41 * mag(u[1].component(2)-u.mean().component(2)),
            dimensionedScalar("", dimensionSet(0,1,-1,0,0), SMALL)
        )
    )
);

massTransferredDwDtPdfTransfer2[0] = sigma[0] * 0.5 *
(
    1 + Foam::erf
    (
        (u[0].component(2) + dt*DuDt0.component(2)) /
        max
        (
            0.5*1.41 * mag(u[0].component(2)-u.mean().component(2)),
            dimensionedScalar("", dimensionSet(0,1,-1,0,0), SMALL)
        )
    )
);
massTransferredDwDtPdfTransfer2[1] = sigma[1] * 0.5 *
(
    1 - Foam::erf
    (
        (u[1].component(2) + dt*DuDt1.component(2)) /
        max
        (
            1*1.41 * mag(u[1].component(2)-u.mean().component(2)),
            dimensionedScalar("", dimensionSet(0,1,-1,0,0), SMALL)
        )
    )
);





massTransferredHorizontalDiv[0] -= dt*max(sign(u[0].component(vector::Z)), dimlessZero)*
min
(
    fvc::div(sigmaf[0]*volFlux[0] - sigmaf[0]*volFlux[0]*mesh.Sf().component(2)/mesh.magSf()),
    rateZero
);
massTransferredHorizontalDiv[1] -= dt*min(sign(u[1].component(vector::Z)), dimlessZero)*
max
(
    fvc::div(sigmaf[1]*volFlux[1] - sigmaf[1]*volFlux[1]*mesh.Sf().component(2)/mesh.magSf()),
    rateZero
);

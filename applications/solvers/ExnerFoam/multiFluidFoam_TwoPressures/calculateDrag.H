// Calculating the drag force between two partitions

if (mag(dragCoeff) > SMALL)
{
    // This code only works for nParts == 2
    if (nParts != 2)
    {
        FatalErrorIn("calculateDrag.H")
            << " drag is only formulated for nParts = 2, not nParts = "
            << nParts << exit(FatalError);
    }

    // Calculate the cloud radius and hence the drag
    cloudRadius = max(cloudRadiusMin, cloudRadiusMax*sigmaf[0]*sigmaf[1]);
    dragCommon = dragCoeff*fvc::interpolate(rho.sum(), "sigmaRho")/cloudRadius
             * mag(Uf[0] - Uf[1]);
    
    // Limit drag for stability
    dragCommon = min
    (
        dragCommon,
        min
        (
            .25/dt*fvc::interpolate(rho[1], "sigmaRho")/sigmaf[0],
            .25/dt*fvc::interpolate(rho[0], "sigmaRho")/sigmaf[1]
        )
    );
}


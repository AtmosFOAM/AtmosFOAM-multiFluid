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
    cloudRadius = max(cloudRadiusMin, cloudRadiusMax*sigma[0]*sigma[1]);
    dragCommon = dragCoeff/cloudRadius * mag(u[0] - u[1]);
}

// Add the heat and momentum transfers
for(label ip = 0; ip < nParts; ip++)
{
    transfer[ip] += dt*dragCommon*sigma[ip];
}
transfer.updateSum();

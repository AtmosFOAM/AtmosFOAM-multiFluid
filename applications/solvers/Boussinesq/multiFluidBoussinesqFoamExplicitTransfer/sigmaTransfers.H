// Mass transfers
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            sigma[ip].oldTime() += dt*
            (
                massTransfer[jp][ip] - massTransfer[ip][jp]
            );
        }
    }
}

// Calculate sigma on the face
for(label ip = 0; ip < nParts; ip++)
{
    sigmaf[ip].oldTime() = fvc::interpolate(sigma[ip].oldTime(), "sigma");
}


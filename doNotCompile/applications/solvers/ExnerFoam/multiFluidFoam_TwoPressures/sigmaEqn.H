// update sigma according to the sigma equation
for(label ip = 0; ip < nParts; ip++)
{
    sigma.ddt()[ip] = -u[ip] & fvc::grad(sigma[ip]);
    sigma[ip] = sigma[ip].oldTime() + dt*sigma.ddt()[ip];
    sigmaf[ip] = fvc::interpolate(sigma[ip], "sigmaRho");
}

sigma.updateSum();
sigmaf.updateSum();


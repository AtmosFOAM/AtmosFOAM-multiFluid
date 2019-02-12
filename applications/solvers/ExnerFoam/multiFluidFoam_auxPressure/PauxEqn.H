// Solve the transport equation for artificial compressibility for each fluid
for(label ip = 0; ip < nParts; ip++)
{
    Paux.ddt()[ip] = -fvc::div(volFlux[ip], Paux[ip], "divPaux");
    Paux[ip] = Paux[ip].oldTime() + dt*
    (
        (1-offCentre)*Paux.ddt()[ip].oldTime()
      + offCentre*Paux.ddt()[ip]
    );
}
Paux.updateSum();
Info << "Paux[0] goes from " << min(Paux[0].internalField()).value()
     << " to " << max(Paux[0].internalField()).value() << endl;
Info << "Paux[1] goes from " << min(Paux[1].internalField()).value()
     << " to " << max(Paux[1].internalField()).value() << endl;
Info << "Paux mean goes from " << min(Paux.mean().internalField()).value()
     << " to " << max(Paux.mean().internalField()).value() << endl;


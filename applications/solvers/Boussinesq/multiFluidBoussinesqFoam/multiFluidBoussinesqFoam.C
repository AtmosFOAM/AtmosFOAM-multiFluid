/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    multiFluidBoussinesqFoam

Description
    Transient Solver for dry, multi-fluid Boussinesq equations

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "PartitionedFields.H"
#include "fvcCurlf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "zeros.H"
    #include "readTransferCoeffs.H"
    #include "readEnvironmentalProps.H"
    #define dt runTime.deltaT()
    #include "createFields.H"
    
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nOuterCorr = itsDict.lookupOrDefault<int>("nOuterCorrectors", 2);
    const int nCorr = itsDict.lookupOrDefault<int>("nCorrectors", 1);
    const int nNonOrthCorr =
        itsDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
    const scalar offCentre = readScalar(mesh.schemesDict().lookup("offCentre"));

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    Info << "sigma[1] goes from " << min(sigma[1].internalField()).value()
         << " to "  << max(sigma[1].internalField()).value() << endl;
    Info << "sigmaf.sum goes from " << min(sigmaf.sum()).value()
         << " to "  << max(sigmaf.sum()).value() << endl;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "partitionedCourantNo.H"

        for (int ucorr=0; ucorr < nOuterCorr; ucorr++)
        {
            #include "sigmaEqn.H"
            if (!noTransfers)
            {
                #include "massTransfers.H"
                #include "applyMassTransfer.H"
            }
            #include "calculateDrag.H"
            #include "bEqn.H"
            #include "PiEqn.H"
//            #include "PEqn.H"
            #include "momentumTransfers.H"
        }

        // Update diagnositcs
        for(label ip = 0; ip < nParts; ip++)
        {
            //sigmaf[ip] = linearInterpolate(sigma[ip]);
            u[ip] = fvc::reconstruct(volFlux[ip]);
            Uf[ip] = linearInterpolate(u[ip]);
            Uf[ip] += (volFlux[ip] - (Uf[ip] & mesh.Sf()))
                      *mesh.Sf()/sqr(mesh.magSf());
            divu[ip] = fvc::div(sigmaf[ip]*volFlux[ip]);
        }
        Uf.updateSum();
        sigmaf.updateSum();
        u.updateSum();
        volFlux.updateSum();
        volVectorField uSum("uSum", fvc::reconstruct(volFlux.mean()));
        uSum.write();
        divu.updateSum();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

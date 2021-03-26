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
    multiFluidBoussinesqFoamPressure.
    A new approach for calculating the pressure in each fluid and the mass
    transfers. As well as the usual multi-fluid equations we solve:
    sigma div(u) = sum(massTransers)
    M_ij = 0.5(simga_i max(-divu_i,0) + sigma_j max(divu_j,0))

Description
    Transient Solver for dry, multi-fluid Boussinesq equations

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "PartitionedFields.H"
#include "TransferFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "zeros.H"
    #include "readEnvironmentalProps.H"
    #include "readTransferCoeffs.H"
    #define dt runTime.deltaT()
    #include "createFields.H"
    
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nOuterCorr = itsDict.lookupOrDefault<int>("nOuterCorrectors", 2);
    const int nNonOrthCorr =
        itsDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
    scalar offCentre = readScalar(mesh.schemesDict().lookup("offCentre"));

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // offCentre = 1 for the first time step
    const scalar offCentreSave = offCentre;
    offCentre = 1;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "partitionedCourantNo.H"

        for (int ucorr=0; ucorr < nOuterCorr; ucorr++)
        {
            if (nParts > 1)
            {
                #include "sigmaEqn.H"
            }
            #include "bEqn.H"
            // Momentum equation with explicit pressure gradient
            #include "momentumEqn.H"

            // Mass transfers
            if (transferType != noTransfer && nParts > 1)
            {
                #include "diffusionTransfers.H"

                #include "massTransfers.H"
                sigma.transferMass(massTransfer, dt);
                interpolate(sigmaf, sigma);
            }
            // Pressure in each fluid
            #include "PEqn.H"
            if (transferType != noTransfer && nParts > 1)
            {
                #include "bTransfers.H"
                #include "momentumTransfers.H"
            }
        }

        Info << "sigma[0] goes from " << min(sigma[0]).value() <<  " to "
            << max(sigma[0]).value() << endl;

        runTime.write();
        offCentre = offCentreSave;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

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
    multiFluidHeatTransferFoam

Description
    Calculate the heat transfer from all fluids and in all directions from
    advection and diffusion

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "PartitionedFields.H"
#include "wallDist.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::addOption
    (
        "region",
        "meshRegion",
        "specify a non-default region to plot"
    );
    argList::addOption
    (
        "cellSet", "cellSetName", "only calculate sums for a subset of cells"
    );
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMeshRegion.H"
    const surfaceScalarField kdir = mesh.Sf().component(2)/mesh.magSf();
    const surfaceScalarField magk = mag(kdir);
    #include "../readEnvironmentalProps.H"
    #include "../readTransferCoeffs.H"
    #define dt runTime.deltaT()
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();
        #include "createFields.H"

        for(label ip = 0; ip < nParts; ip++)
        {
            heatTransferf[ip] = volFlux[ip]*linearInterpolate(b[ip])/mesh.magSf()
                              - alphaf*fvc::snGrad(b[ip]);
        }
        if (maxAlpha0.value() > SMALL)
        {
            heatTransferf[0] -= alpha0*fvc::snGrad(b[0]);
        }
        for(label ip = 0; ip < nParts; ip++)
        {
            heatTransfer[ip] = fvc::reconstruct(heatTransferf[ip]*mesh.magSf());
/*            forAll(heatTransfer[ip].boundaryField(), ipat)
            {
                heatTransfer[ip].boundaryFieldRef()[ipat]
                     = heatTransferf[ip].boundaryField()[ipat]
                      *mesh.boundary()[ipat].Sf()/mesh.boundary()[ipat].magSf();
            }
*/        }
        heatTransfer.write();
        heatTransferf.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

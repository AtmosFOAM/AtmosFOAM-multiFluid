/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

InClass
    Foam::TransferField

Description

SourceFiles
    TransferFields.C

\*---------------------------------------------------------------------------*/

#ifndef TransferFieldsFwd_H
#define TransferFieldsFwd_H

#include "TransferField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class volMesh;
class surfaceMesh;

template<class Type>
class fvPatchField;

template<class Type>
class fvsPatchField;

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricField;

typedef TransferField<scalar, fvPatchField, volMesh> volScalarTransferField;
typedef TransferField<scalar, fvsPatchField, surfaceMesh> surfaceScalarTransferField;
typedef TransferField<vector, fvPatchField, volMesh> volVectorTransferField;
typedef TransferField<vector, fvsPatchField, surfaceMesh> surfaceVectorTransferField;

template<class Type>
inline TransferField<Type, fvsPatchField, surfaceMesh>& interpolate
(
    TransferField<Type, fvsPatchField, surfaceMesh>& Tf,
    const TransferField<Type, fvPatchField, volMesh>& T
)
{
    for(label ip = 0; ip < Tf.size(); ip++)
    {
        Tf[ip] = linearInterpolate(T[ip]);
    }
    return Tf;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::PartitionedField

Description

SourceFiles
    PartitionedFields.C

\*---------------------------------------------------------------------------*/

#ifndef PartitionedFieldsFwd_H
#define PartitionedFieldsFwd_H

#include "PartitionedField.H"

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

typedef PartitionedField<scalar, fvPatchField, volMesh> partitionedVolScalarField;
typedef PartitionedField<scalar, fvsPatchField, surfaceMesh> partitionedSurfaceScalarField;
typedef PartitionedField<vector, fvPatchField, volMesh> partitionedVolVectorField;
typedef PartitionedField<vector, fvsPatchField, surfaceMesh> partitionedSurfaceVectorField;

template<class Type>
inline PartitionedField<Type, fvsPatchField, surfaceMesh>& interpolate
(
    PartitionedField<Type, fvsPatchField, surfaceMesh>& Tf,
    const PartitionedField<Type, fvPatchField, volMesh>& T
)
{
    for(label ip = 0; ip < Tf.size(); ip++)
    {
        Tf[ip] = linearInterpolate(T[ip]);
    }
    Tf.updateSum();
    return Tf;
}

template<class Type>
inline PtrList<PartitionedField<Type, fvsPatchField, surfaceMesh>>& interpolate
(
    PtrList<PartitionedField<Type, fvsPatchField, surfaceMesh>>& Tf,
    const PtrList<PartitionedField<Type, fvPatchField, volMesh>>& T
)
{
    for(label ip = 0; ip < Tf.size(); ip++)
    {
        for(label jp = 0; jp < Tf[ip].size(); jp++)
        {
            Tf[ip][jp] = linearInterpolate(T[ip][jp]);
        }
        Tf[ip].updateSum();
    }
    return Tf;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

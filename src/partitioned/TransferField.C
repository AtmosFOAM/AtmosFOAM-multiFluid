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

\*---------------------------------------------------------------------------*/

#include "TransferField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::TransferField<Type, PatchField, GeoMesh>::TransferField
(
    const word& baseName__,
    const wordList& partNames__,
    const Mesh& mesh,
    const word& timeName,
    IOobject::writeOption writeOpt
)
:
    PtrList<GeometricField<Type, PatchField, GeoMesh> >
    (
        partNames__.size()*(partNames__.size()-1)
    ),
    baseName_(baseName__),
    partNames_(partNames__)
{
    label it = 0;
    for(label ip = 0; ip < nParts(); ip++)
    {
        for(label jp = 0; jp < nParts(); jp++)
        {
            if (ip != jp)
            {
                this->set
                (
                    it,
                    new GeometricField<Type, PatchField, GeoMesh>
                    (
                        IOobject
                        (
                            baseName()+'.'+partNames()[ip]+'.'+partNames()[jp],
                            timeName,
                            mesh,
                            IOobject::MUST_READ,
                            writeOpt
                        ),
                        mesh
                    )
                );
//                Info << "ip = " << ip << " jp = " << jp << " it = " << it
//                     << endl;
                it++;
            }
        }
    }
}

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::TransferField<Type, PatchField, GeoMesh>::TransferField
(
    const word& baseName__,
    const wordList& partNames__,
    const GeometricField<Type, PatchField, GeoMesh>& field,
    IOobject::writeOption writeOpt
)
:
    PtrList<GeometricField<Type, PatchField, GeoMesh> >
    (
        partNames__.size()*(partNames__.size()-1)
    ),
    baseName_(baseName__),
    partNames_(partNames__)
{
    label it = 0;
    for(label ip = 0; ip < nParts(); ip++)
    {
        for(label jp = 0; jp < nParts(); jp++)
        {
            if (ip != jp)
            {
                this->set
                (
                    it,
                    new GeometricField<Type, PatchField, GeoMesh>
                    (
                        IOobject
                        (
                            baseName()+'.'+partNames()[ip]+'.'+partNames()[jp],
                           field.mesh().time().timeName(),
                           field.mesh(),
                           IOobject::READ_IF_PRESENT,
                           writeOpt
                       ),
                       field
                   )
                );
//                Info << "ip = " << ip << " jp = " << jp << " it = " << it
//                     << endl;
                it++;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::TransferField<Type, PatchField, GeoMesh>::~TransferField()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::TransferField<Type, PatchField, GeoMesh>::write()
{
    for(label ip = 0; ip < this->size(); ip++)
    {
        this->operator[](ip).write();
    }

}

// * * * * * * * * * * * * * Functions from PtrList * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
const Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::TransferField<Type, PatchField, GeoMesh>::operator()
(
    const label ip, label jp
) const
{
    if (ip == jp || ip < 0 || jp < 0 || ip >= nParts() || jp >= nParts())
    {
        FatalErrorIn("TransferField")
             << " cannot access transfer field for ip = " << ip
             << " jp = " << jp << exit(FatalError);
    }
    label j = jp < ip ? jp : (jp > ip ? jp-1 : -1);
    label it = ip*(nParts()-1) + j;
    return PtrList<GeometricField<Type, PatchField, GeoMesh> >::operator[](it);
}

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::TransferField<Type, PatchField, GeoMesh>::operator()
(
    const label ip, label jp
)
{
    if (ip == jp || ip < 0 || jp < 0 || ip >= nParts() || jp >= nParts())
    {
        FatalErrorIn("TransferField")
             << " cannot access transfer field for ip = " << ip
             << " jp = " << jp << exit(FatalError);
    }
    label j = jp < ip ? jp : (jp > ip ? jp-1 : -1);
    label it = ip*(nParts()-1) + j;
//    Info << "ip = " << ip << " jp = " << jp << " it = " << it << endl;
    return PtrList<GeometricField<Type, PatchField, GeoMesh> >::operator[](it);
}

// ************************************************************************* //

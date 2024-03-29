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

#include "PartitionedField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::PartitionedField
(
    const word& baseName__,
    const wordList& partNames__,
    const Mesh& mesh,
    const word& timeName,
    IOobject::writeOption writeOpt
)
:
    PtrList<GeometricField<Type, PatchField, GeoMesh> >(partNames__.size()),
    baseName_(baseName__),
    partNames_(partNames__),
    sum_
    (
        IOobject
        (
            baseName_+".sum", timeName, mesh,
            IOobject::NO_READ, writeOpt
        ),
        mesh,
        dimensioned<Type>("sum", dimless, pTraits<Type>::one)
    ),
    mean_
    (
        IOobject
        (
            baseName_, timeName, mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<Type>("mean", dimless, pTraits<Type>::one)
    ),
    needsSigma_(false),
    sigmaPtr_(NULL),
    ddtPtr_(NULL)
{
    for(label ip = 0; ip < size(); ip++)
    {
        (*this).set
        (
            ip,
            new GeometricField<Type, PatchField, GeoMesh>
            (
                IOobject
                (
                    baseName() + '.' + partNames()[ip], timeName, mesh,
                    IOobject::MUST_READ, IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }

    updateMean();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::PartitionedField
(
    const word& baseName__,
    const wordList& partNames__,
    const Mesh& mesh,
    const word& timeName,
    const PartitionedField<scalar, PatchField, GeoMesh>& sigma__,
    IOobject::writeOption writeOpt
)
:
    PtrList<GeometricField<Type, PatchField, GeoMesh> >(partNames__.size()),
    baseName_(baseName__),
    partNames_(partNames__),
    sum_
    (
        IOobject
        (
            baseName_+".sum", timeName, mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<Type>("sum", dimless, pTraits<Type>::one)
    ),
    mean_
    (
        IOobject
        (
            baseName_, timeName, mesh,
            IOobject::NO_READ, writeOpt
        ),
        mesh,
        dimensioned<Type>("mean", dimless, pTraits<Type>::one)
    ),
    needsSigma_(true),
    sigmaPtr_(&sigma__),
    ddtPtr_(NULL)
{
    for(label ip = 0; ip < size(); ip++)
    {
        (*this).set
        (
            ip,
            new GeometricField<Type, PatchField, GeoMesh>
            (
                IOobject
                (
                    baseName() + '.' + partNames()[ip], timeName, mesh,
                    IOobject::MUST_READ, IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }

    updateMean();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::PartitionedField
(
    const word& baseName__,
    const wordList& partNames__,
    const GeometricField<Type, PatchField, GeoMesh>& field,
    IOobject::writeOption writeOpt
)
:
    PtrList<GeometricField<Type, PatchField, GeoMesh> >(partNames__.size()),
    baseName_(baseName__),
    partNames_(partNames__),
    sum_
    (
        IOobject
        (
            baseName_+".sum", field.time().timeName(), field.mesh(),
            IOobject::NO_READ, writeOpt
        ),
        field.mesh(),
        dimensioned<Type>("sum", dimless, pTraits<Type>::one)
    ),
    mean_
    (
        IOobject
        (
            baseName_, field.time().timeName(), field.mesh(),
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        field.mesh(),
        dimensioned<Type>("mean", dimless, pTraits<Type>::one)
    ),
    needsSigma_(false),
    sigmaPtr_(NULL),
    ddtPtr_(NULL)
{
    for(label ip = 0; ip < size(); ip++)
    {
        (*this).set
        (
            ip,
            new GeometricField<Type, PatchField, GeoMesh>
            (
                IOobject
                (
                    baseName() + '.' + partNames()[ip],
                    field.mesh().time().timeName(),
                    field.mesh(),
                    IOobject::READ_IF_PRESENT, writeOpt
                ),
                field
            )
        );
    }
    updateMean();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::PartitionedField
(
    const word& baseName__,
    const wordList& partNames__,
    const GeometricField<Type, PatchField, GeoMesh>& field,
    const PartitionedField<scalar, PatchField, GeoMesh>& sigma__,
    IOobject::writeOption writeOpt
)
:
    PtrList<GeometricField<Type, PatchField, GeoMesh> >(partNames__.size()),
    baseName_(baseName__),
    partNames_(partNames__),
    sum_
    (
        IOobject
        (
            baseName_+".sum", field.time().timeName(), field.mesh(),
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        field.mesh(),
        dimensioned<Type>("sum", dimless, pTraits<Type>::one)
    ),
    mean_
    (
        IOobject
        (
            baseName_, field.time().timeName(), field.mesh(),
            IOobject::NO_READ, writeOpt
        ),
        field.mesh(),
        dimensioned<Type>("mean", dimless, pTraits<Type>::one)
    ),
    needsSigma_(true),
    sigmaPtr_(&sigma__),
    ddtPtr_(NULL)
{
    for(label ip = 0; ip < size(); ip++)
    {
        (*this).set
        (
            ip,
            new GeometricField<Type, PatchField, GeoMesh>
            (
                IOobject
                (
                    baseName_+'.'+partNames_[ip],
                    field.mesh().time().timeName(),
                    field.mesh(),
                    IOobject::READ_IF_PRESENT, writeOpt
                ),
                field
            )
        );
    }
    
    updateMean();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::PartitionedField
(
    const PartitionedField& pf__
)
:
    PtrList<GeometricField<Type, PatchField, GeoMesh> >(pf__),
    baseName_(pf__.baseName_),
    partNames_(pf__.partNames_),
    sum_(pf__.sum_),
    mean_(pf__.mean_),
    needsSigma_(pf__.needsSigma_),
    sigmaPtr_(pf__.sigmaPtr_),
    ddtPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::~PartitionedField()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
const Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::PartitionedField<Type, PatchField, GeoMesh>::updateSum()
{
    mean_.dimensions().reset(operator[](0).dimensions());

    // Depends on whether it is sigma weighted
    
    if (!needsSigma_)
    {
        sum_.dimensions().reset(operator[](0).dimensions());
        sum_ = operator[](0);

        // Sum contributions from other partitions
        for (label ip = 1; ip < size(); ip++)
        {
            sum_ += operator[](ip);
        }
        mean_ = sum_;
    }
    else
    {
        sum_.dimensions().reset
        (
            operator[](0).dimensions()*sigma()[0].dimensions()
        );
    
        sum_ = sigma()[0]*operator[](0);

        // Sum contributions from other partitions
        for (label ip = 1; ip < size(); ip++)
        {
            sum_ += sigma()[ip]*operator[](ip);
        }
        mean_ = sum_/
            max
            (
                sigma().sum(),
                dimensionedScalar("", sigma().sum().dimensions(), SMALL)
            );
    }
    
    return sum_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
const Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::PartitionedField<Type, PatchField, GeoMesh>::updateMean()
{
    updateSum();
    return mean_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::timesSigma() const
{
    if (!needsSigma())
    {
        FatalErrorIn("PartitionedField::timesSigma")
            << "cannot calculate timesSigma without sigma"
            << abort(FatalError);
    }
    
    // Name of the timesSigma
    string newName = baseName()+'.'+sigma().baseName();
    
    // Geometric field to start from
    const GeometricField<Type, PatchField, GeoMesh>& f = operator[](0);
    GeometricField<Type, PatchField, GeoMesh> fracField
    (
        IOobject(newName, f.time().timeName(), f.mesh()),
        f*sigma()[0],
        f.boundaryField().types()
    );
    
    PartitionedField<Type, PatchField, GeoMesh> frac
    (
        newName, partNames(), fracField
    );
    
    for(label ip = 1; ip < size(); ip++)
    {
        frac[ip] = operator[](ip)*sigma()[ip];
    }
    
    frac.updateSum();
    
    return frac;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::divideBy
(
    const PartitionedField<scalar, PatchField, GeoMesh>& sigma__,
    const word newBaseName
) const
{
    if (needsSigma())
    {
        FatalErrorIn("PartitionedField::divideBy")
            << "cannot calculate divideBy sigma for a field with sigma"
            << abort(FatalError);
    }

    // Name of the new field
    string newName = newBaseName;
    if (newBaseName == "")
    {
        newName = baseName();
        newName.erase(sigma()[0].name().size(), newName.size()-1);
    }
    
    // Geometric field to start from
    GeometricField<Type, PatchField, GeoMesh> ff
    (
        newName, operator[](0)/sigma__[0]
    );

    // New partitionedField
    PartitionedField<Type, PatchField, GeoMesh> f
    (
        newName, partNames(), ff, sigma__
    );
    
    for(label ip = 1; ip < size(); ip++)
    {
        f[ip] = operator[](ip)/sigma__[ip];
    }

    return f;
}

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::PartitionedField<Type, PatchField, GeoMesh>::transferMass
(
    const TransferField<Type, PatchField, GeoMesh>& M,
    const dimensionedScalar& dt
)
{
    for(label ip = 0; ip < size(); ip++)
    {
        for(label jp=0; jp < size(); jp++) if (ip != jp)
        {
            operator[](ip) += dt*(M(jp,ip) - M(ip,jp));
        }
        operator[](ip).correctBoundaryConditions();
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::PartitionedField<Type, PatchField, GeoMesh>::transferField
(
    const TransferField<Type, PatchField, GeoMesh>& M,
    const TransferField<Type, PatchField, GeoMesh>& fieldT,
    const dimensionedScalar& dt
)
{
    const dimensionedScalar smallSigma("", sigma()[0].dimensions(), SMALL);
    const PartitionedField<Type, PatchField, GeoMesh> old = *this;

    for(label ip = 0; ip < size(); ip++)
    {
        for(label jp = 0; jp < size(); jp++) if (ip != jp)
        {
            operator[](ip) += dt/max(sigma()[ip], smallSigma)*
            (
                M(jp,ip)*(fieldT(jp,ip) - old[ip])
              - M(ip,jp)*(fieldT(ip,jp) - old[ip])
            );
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::PartitionedField<Type, PatchField, GeoMesh>::transferField
(
    const TransferField<Type, PatchField, GeoMesh>& M,
    const dimensionedScalar& dt
)
{
    const dimensionedScalar smallSigma("", sigma()[0].dimensions(), SMALL);
    const PartitionedField<Type, PatchField, GeoMesh> old = *this;

    for(label ip = 0; ip < size(); ip++)
    {
        for(label jp = ip+1; jp < size(); jp++)
        {
            operator[](ip) += dt/max(sigma()[ip], smallSigma)*
            (
                M(jp,ip)*(old[jp] - old[ip])
            );
            operator[](jp) += dt/max(sigma()[jp], smallSigma)*
            (
                M(ip,jp)*(old[ip] - old[jp])
            );
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::PartitionedField<Type, PatchField, GeoMesh>
         ::asymetricTransferField
(
    const TransferField<Type, PatchField, GeoMesh>& M,
    const TransferField<Type, PatchField, GeoMesh>& fieldT,
    const TransferField<Type, PatchField, GeoMesh>& add,
    const TransferField<Type, PatchField, GeoMesh>& remove,
    const dimensionedScalar& dt
)
{
    const dimensionedScalar smallSigma("", sigma()[0].dimensions(), SMALL);
    const PartitionedField<Type, PatchField, GeoMesh> old = *this;

    for(label ip = 0; ip < size(); ip++)
    {
        for(label jp = 0; jp < size(); jp++) if (ip != jp)
        {
            operator[](ip) += dt/max(sigma()[ip], smallSigma)*
            (
                M(jp,ip)*(fieldT(jp,ip) - old[ip] + add(jp,ip))
              - M(ip,jp)*(fieldT(ip,jp) - old[ip] + remove(ip,jp))
            );
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::PartitionedField<Type, PatchField, GeoMesh>::storeTime()
{
    // Only call this function once (only set the pointers once)
    if (ddtPtr_)
    {
        FatalErrorIn("PartitionedField::storeTime")
            << "already called, ddtPtr_ already set"
            << abort(FatalError);
    }

    // Store old time for all of the fields
    for(label ip = 0; ip < size(); ip++)
    {
        operator[](ip).oldTime();
    }
    
    // Store old sum and mean
    sum().oldTime();
    mean().oldTime();
    
    // Create fields for all the time directories
    ddtPtr_ = new PartitionedField<Type, PatchField, GeoMesh>
    (
        "ddt"+baseName(),
        partNames(),
        GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                "ddt."+baseName()+"."+partNames()[0],
                operator[](0).mesh().time().timeName(),
                operator[](0).mesh()
            ),
            (operator[](0) - operator[](0).oldTime())
                /operator[](0).time().deltaT(),
            "fixedValue"
        ),
        IOobject::NO_WRITE
    );
    ddt()[0].oldTime();
    for(label ip = 1; ip < size(); ip++)
    {
        ddt()[ip] = (operator[](ip) - operator[](ip).oldTime())
                    /operator[](ip).time().deltaT();
        ddt()[ip].oldTime();
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::PartitionedField<Type, PatchField, GeoMesh>::write()
{
    for(label ip = 0; ip < size(); ip++)
    {
        operator[](ip).write();
    }

    updateSum();
    mean_.write();

    // Write out old time and rate of change if set
    if (ddtPtr_)
    {
        ddt().write();
        for(label ip = 0; ip < size(); ip++)
        {
            operator[](ip).oldTime().write();
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::PartitionedField<Type, PatchField, GeoMesh>::readUpdate()
{
    for(label ip = 0; ip < size(); ip++)
    {
        operator[](ip).readUpdate();
    }
    updateSum();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::PartitionedField<Type, PatchField, GeoMesh>::operator=
(
    const PartitionedField<Type, PatchField, GeoMesh>& gf
)
{
    if (this == &gf)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }
    
    for(label ip = 0; ip < size(); ip++)
    {
        operator[](ip) = gf[ip];
    }
    sum_ = gf.sum();
    mean_ = gf.mean();
}

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::PartitionedField<Type, PatchField, GeoMesh>::operator+=
(
    const PartitionedField<Type, PatchField, GeoMesh>& gf
)
{
    // Check that you are adding two partitions with the same sigma
    if (sigmaPtr_ != &gf.sigma())
    {
        FatalErrorIn("PartitionedField::operator+=")
            << " attempting to add two fields with different sigmas"
            << abort(FatalError);
    }

    for(label ip = 0; ip < size(); ip++)
    {
        operator[](ip) += gf[ip];
    }
    
    updateSum();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::PartitionedField<Type, PatchField, GeoMesh>::operator-=
(
    const PartitionedField<Type, PatchField, GeoMesh>& gf
)
{
    // Check that you are subtracting two partitions with the same sigma
    if (sigmaPtr_ != &gf.sigma())
    {
        FatalErrorIn("PartitionedField::operator-=")
            << " attempting to add two fields with different sigmas"
            << abort(FatalError);
    }

    for(label ip = 0; ip < size(); ip++)
    {
        operator[](ip) -= gf[ip];
    }

    updateSum();
}


// ************************************************************************* //

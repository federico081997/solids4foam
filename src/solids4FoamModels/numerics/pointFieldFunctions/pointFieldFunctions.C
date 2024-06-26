/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "pointFieldFunctions.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

tmp<pointSymmTensorField> symm(const pointTensorField& ptr)
{
	const pointMesh& pMesh =  ptr.mesh();

	// Prepare the temporary field
	tmp<pointSymmTensorField> tresult
	(
	    new pointSymmTensorField
	    (
	        IOobject
	        (
	            "symmetricField",
	            pMesh.mesh().time().timeName(),
	            pMesh.mesh(),
	            IOobject::NO_READ,
	            IOobject::NO_WRITE
	        ),
	        pMesh,
	        dimensionedSymmTensor("0", dimless, symmTensor::zero)
	    )
	);

#ifdef OPENFOAM_NOT_EXTEND
	pointSymmTensorField& result = tresult.ref();
#else
	pointSymmTensorField& result = tresult();
#endif

	// Set the result field to be symm(ptr)
#ifdef OPENFOAM_NOT_EXTEND
    result.primitiveFieldRef() = symm(ptr.primitiveField());
#else
    result.internalField() = symm(ptr.internalField());
#endif

	result.correctBoundaryConditions();

	return tresult;
}

tmp<pointSymmTensorField> dev(const pointSymmTensorField& ptr)
{
	const pointMesh& pMesh =  ptr.mesh();

	// Prepare the temporary field
	tmp<pointSymmTensorField> tresult
	(
	    new pointSymmTensorField
	    (
	        IOobject
	        (
	            "devField",
	            pMesh.mesh().time().timeName(),
	            pMesh.mesh(),
	            IOobject::NO_READ,
	            IOobject::NO_WRITE
	        ),
	        pMesh,
	        dimensionedSymmTensor("0", dimless, symmTensor::zero)
	    )
	);

#ifdef OPENFOAM_NOT_EXTEND
	pointSymmTensorField& result = tresult.ref();
#else
	pointSymmTensorField& result = tresult();
#endif

	// Set the result field to be dev(ptr)
#ifdef OPENFOAM_NOT_EXTEND
    result.primitiveFieldRef() = dev(ptr.primitiveField());
#else
    pointD.internalField() = dev(ptr.internalField());
#endif

	result.correctBoundaryConditions();

	return tresult;
}

tmp<pointScalarField> tr(const pointSymmTensorField& ptr)
{
	const pointMesh& pMesh =  ptr.mesh();

	// Prepare the temporary field
	tmp<pointScalarField> tresult
	(
	    new pointScalarField
	    (
	        IOobject
	        (
	            "trField",
	            pMesh.mesh().time().timeName(),
	            pMesh.mesh(),
	            IOobject::NO_READ,
	            IOobject::NO_WRITE
	        ),
	        pMesh,
	        dimensionedScalar("0", dimless, 0.0)
	    )
	);

#ifdef OPENFOAM_NOT_EXTEND
	pointScalarField& result = tresult.ref();
#else
	pointScalarField& result = tresult();
#endif

	// Set the result field to be tr(ptr)
#ifdef OPENFOAM_NOT_EXTEND
    result.primitiveFieldRef() = tr(ptr.primitiveField());
#else
    pointD.internalField() = tr(ptr.internalField());
#endif

	result.correctBoundaryConditions();

	return tresult;
}

tmp<pointTensorField> cof(const pointTensorField& ptr)
{
	const pointMesh& pMesh =  ptr.mesh();

	// Prepare the temporary field
	tmp<pointTensorField> tresult
	(
	    new pointTensorField
	    (
	        IOobject
	        (
	            "cofField",
	            pMesh.mesh().time().timeName(),
	            pMesh.mesh(),
	            IOobject::NO_READ,
	            IOobject::NO_WRITE
	        ),
	        pMesh,
	        dimensionedTensor("zero", dimless, tensor::zero)
	    )
	);

#ifdef OPENFOAM_NOT_EXTEND
	pointTensorField& result = tresult.ref();
#else
	pointTensorField& result = tresult();
#endif

	// Set the result field to be cof(ptr)
#ifdef OPENFOAM_NOT_EXTEND
    result.primitiveFieldRef() = cof(ptr.primitiveField());
#else
    pointD.internalField() = cof(ptr.internalField());
#endif

	result.correctBoundaryConditions();

	return tresult;
}

tmp<pointScalarField> sqr(const pointScalarField& ptr)
{
	const pointMesh& pMesh =  ptr.mesh();

	// Prepare the temporary field
	tmp<pointScalarField> tresult
	(
	    new pointScalarField
	    (
	        IOobject
	        (
	            "sqrField",
	            pMesh.mesh().time().timeName(),
	            pMesh.mesh(),
	            IOobject::NO_READ,
	            IOobject::NO_WRITE
	        ),
	        pMesh,
	        dimensionedScalar("0", dimless, 0.0)
	    )
	);

#ifdef OPENFOAM_NOT_EXTEND
	pointScalarField& result = tresult.ref();
#else
	pointScalarField& result = tresult();
#endif

	// Set the result field to be sqr(ptr)
#ifdef OPENFOAM_NOT_EXTEND
    result.primitiveFieldRef() = sqr(ptr.primitiveField());
#else
    pointD.internalField() = sqr(ptr.internalField());
#endif

	result.correctBoundaryConditions();

	return tresult;
}

tmp<pointScalarField> pow
(
    const pointScalarField& ptr,
    const scalar& power
)
{
	const pointMesh& pMesh =  ptr.mesh();

	// Prepare the temporary field
	tmp<pointScalarField> tresult
	(
	    new pointScalarField
	    (
	        IOobject
	        (
	            "powField",
	            pMesh.mesh().time().timeName(),
	            pMesh.mesh(),
	            IOobject::NO_READ,
	            IOobject::NO_WRITE
	        ),
	        pMesh,
	        dimensionedScalar("0", dimless, 0.0)
	    )
	);

#ifdef OPENFOAM_NOT_EXTEND
	pointScalarField& result = tresult.ref();
#else
	pointScalarField& result = tresult();
#endif

	// Set the result field to be pow(ptr, power)
#ifdef OPENFOAM_NOT_EXTEND
    result.primitiveFieldRef() = pow(ptr.primitiveField(), power);
#else
    pointD.internalField() = pow(ptr.internalField(), power);
#endif

	result.correctBoundaryConditions();

	return tresult;
}

}

// ************************************************************************* //

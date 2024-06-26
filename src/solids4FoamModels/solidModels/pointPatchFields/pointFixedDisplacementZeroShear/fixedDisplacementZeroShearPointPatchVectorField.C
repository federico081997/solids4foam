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

#include "fixedDisplacementZeroShearPointPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "pointPatchFields.H"
#include "pointBoundaryMesh.H"
#include "pointMesh.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "Time.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedDisplacementZeroShearPointPatchVectorField::
fixedDisplacementZeroShearPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    componentMixedPointPatchVectorField(p, iF),
    dispSeries_(),
    curTimeIndex_(-1)
{}


fixedDisplacementZeroShearPointPatchVectorField::
fixedDisplacementZeroShearPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    componentMixedPointPatchVectorField(p, iF),
    dispSeries_(),
    curTimeIndex_(-1)
{
#ifdef OPENFOAM_NOT_EXTEND
    // refValue is only used by foam-extend
    vectorField refValue(patch().size(), vector::zero);
#else
    // Set value fraction to fixed normal direction
    // valueFraction() = sqr(patch().pointNormals());
    valueFraction() = patch().pointNormals();

    vectorField& refValue = this->refValue();
#endif

    // Set displacement
    if (dict.found("dispSeries"))
    {
        if (debug)
        {
            Info<< "    disp is time-varying" << endl;
        }

        dispSeries_ =
            interpolationTable<vector>(dict.subDict("dispSeries"));

        refValue = dispSeries_(db().time().timeOutputValue());
    }
    else
    {
        // Use the "value" on the patch
        refValue = this->patchInternalField();
    }

    //this->updateBoundaryField();

    // Set the boundary values

    tmp<vectorField> internalValues = this->patchInternalField();

#ifdef OPENFOAM_NOT_EXTEND
    this->setInInternalField
    (
        const_cast<vectorField&>(this->primitiveField()),
        refValue
    );
#else
    const vectorField values
    (
        cmptMultiply(refValue, valueFraction())
      + cmptMultiply(internalValues, vector::one - valueFraction())
    );

    this->setInInternalField
    (
        const_cast<vectorField&>(this->internalField()), values
    );
#endif
}


fixedDisplacementZeroShearPointPatchVectorField::
fixedDisplacementZeroShearPointPatchVectorField
(
    const fixedDisplacementZeroShearPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const PointPatchFieldMapper&
)
:
    componentMixedPointPatchVectorField(p, iF),
    dispSeries_(),
    curTimeIndex_(-1)
{
#ifdef FOAMEXTEND
    refValue() = vector::zero;
    valueFraction() = patch().pointNormals();
#endif
}


#ifndef OPENFOAM_ORG
fixedDisplacementZeroShearPointPatchVectorField::
fixedDisplacementZeroShearPointPatchVectorField
(
    const fixedDisplacementZeroShearPointPatchVectorField& ptf
)
:
    componentMixedPointPatchVectorField(ptf),
    dispSeries_(),
    curTimeIndex_(-1)
{}
#endif


fixedDisplacementZeroShearPointPatchVectorField::
fixedDisplacementZeroShearPointPatchVectorField
(
    const fixedDisplacementZeroShearPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    componentMixedPointPatchVectorField(ptf, iF),
    dispSeries_(),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map and resize from self given a mapper
void fixedDisplacementZeroShearPointPatchVectorField::autoMap
(
    const PointPatchFieldMapper& m
)
{
    componentMixedPointPatchVectorField::autoMap(m);
}


// Grab the values using rmap
void fixedDisplacementZeroShearPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    componentMixedPointPatchVectorField::rmap(ptf, addr);
}


void fixedDisplacementZeroShearPointPatchVectorField::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (curTimeIndex_ != db().time().timeIndex())
    {
#ifdef FOAMEXTEND
        // If time-varying, update the displacement field
        if (dispSeries_.size())
        {
            refValue() = dispSeries_(db().time().timeOutputValue());
        }
#endif

        curTimeIndex_ = db().time().timeIndex();
    }

    componentMixedPointPatchVectorField::initEvaluate(commsType);
}


void fixedDisplacementZeroShearPointPatchVectorField::write(Ostream& os) const
{
    componentMixedPointPatchVectorField::write(os);

    if (dispSeries_.size())
    {
        os.writeKeyword("dispSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        dispSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    fixedDisplacementZeroShearPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

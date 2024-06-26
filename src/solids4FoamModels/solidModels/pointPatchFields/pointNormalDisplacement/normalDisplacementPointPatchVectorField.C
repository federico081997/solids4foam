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

#include "normalDisplacementPointPatchVectorField.H"
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

normalDisplacementPointPatchVectorField::
normalDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(p, iF),
    dispSeries_()
{}


normalDisplacementPointPatchVectorField::
normalDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchVectorField(p, iF),
    dispSeries_()
{
    Info<< "Creating " << type() << " boundary condition" << endl;

    // Check if displacement is time-varying
    if (dict.found("displacementSeries"))
    {
        Info<< "    displacement is time-varying" << endl;
        dispSeries_ =
            interpolationTable<scalar>(dict.subDict("displacementSeries"));
    }
}


normalDisplacementPointPatchVectorField::
normalDisplacementPointPatchVectorField
(
    const normalDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const PointPatchFieldMapper&
)
:
    fixedValuePointPatchVectorField(p, iF),
    dispSeries_()
{}


#ifndef OPENFOAM_ORG
normalDisplacementPointPatchVectorField::
normalDisplacementPointPatchVectorField
(
    const normalDisplacementPointPatchVectorField& ptf
)
:
    fixedValuePointPatchVectorField(ptf),
    dispSeries_()
{}
#endif


normalDisplacementPointPatchVectorField::
normalDisplacementPointPatchVectorField
(
    const normalDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(ptf, iF),
    dispSeries_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map and resize from self given a mapper
void normalDisplacementPointPatchVectorField::autoMap
(
    const PointPatchFieldMapper& m
)
{
    fixedValuePointPatchVectorField::autoMap(m);
}


// Grab the values using rmap
void normalDisplacementPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValuePointPatchVectorField::rmap(ptf, addr);
}

void normalDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    this->operator==
    (
        patch().pointNormals()*dispSeries_(this->db().time().timeOutputValue())
    );
    fixedValuePointPatchField<vector>::updateCoeffs();
}


void normalDisplacementPointPatchVectorField::write(Ostream& os) const
{
    fixedValuePointPatchVectorField::write(os);

    os.writeKeyword("displacementSeries") << nl;
    os << token::BEGIN_BLOCK << nl;
    dispSeries_.write(os);
    os << token::END_BLOCK << nl;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    normalDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

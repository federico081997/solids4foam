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

#include "solidTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "lookupSolidModel.H"
#include "patchCorrectionVectors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    nonOrthogonalCorrections_(true),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    tractionSeries_(),
    pressureSeries_(),
    tractionFieldPtr_(),
    pressureFieldPtr_(),
    secondOrder_(false),
    setEffectiveTraction_(false),
    relaxFac_(1.0),
    curTimeIndex_(-1)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    nonOrthogonalCorrections_
    (
        dict.lookupOrDefault<Switch>("nonOrthogonalCorrections", true)
    ),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    tractionSeries_(),
    pressureSeries_(),
    tractionFieldPtr_(),
    pressureFieldPtr_(),
    secondOrder_(dict.lookupOrDefault<Switch>("secondOrder", false)),
    setEffectiveTraction_
    (
        dict.lookupOrDefault<Switch>("setEffectiveTraction", false)
    ),
    relaxFac_(dict.lookupOrDefault<scalar>("relaxationFactor", 1.0)),
    curTimeIndex_(-1)
{
    Info<< "Creating " << type() << " boundary condition" << endl;

    if (dict.found("gradient"))
    {
        gradient() = vectorField("gradient", dict, p.size());
    }
    else
    {
        gradient() = vector::zero;
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }

    // Check how traction is defined
    if
    (
        (dict.found("traction") && dict.found("tractionSeries"))
     || (dict.found("traction") && dict.found("tractionField"))
     || (dict.found("tractionSeries") && dict.found("tractionField"))
    )
    {
        FatalErrorIn
        (
            "solidTractionFvPatchVectorField::solidTractionFvPatchVectorField"
        )   << "Only one of traction, tractionSeries or tractionField can be "
            << "specified!"
            << abort(FatalError);
    }
    else if (dict.found("tractionSeries"))
    {
        Info<< "    traction is time-varying" << endl;
        tractionSeries_ =
            interpolationTable<vector>(dict.subDict("tractionSeries"));
    }
    else if (dict.found("tractionField"))
    {
        Info<< "    traction is specified as a field" << endl;
        tractionFieldPtr_.set
        (
            new volVectorField
            (
                IOobject
                (
                    word(dict.lookup("tractionField")),
                    patch().boundaryMesh().mesh().time().timeName(),
                    patch().boundaryMesh().mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh()
            )
        );
    }
    else
    {
        traction_ = vectorField("traction", dict, p.size());
    }

    // Check how pressure is defined
    if
    (
        (dict.found("pressure") && dict.found("pressureSeries"))
     || (dict.found("pressure") && dict.found("pressureField"))
     || (dict.found("pressureSeries") && dict.found("pressureField"))
    )
    {
        FatalErrorIn
        (
            "solidTractionFvPatchVectorField::solidTractionFvPatchVectorField"
        )   << "Only one of pressure, pressureSeries or pressureField can be "
            << "specified!"
            << abort(FatalError);
    }
    else if (dict.found("pressureSeries"))
    {
        Info<< "    pressure is time-varying" << endl;
        pressureSeries_ =
            interpolationTable<scalar>(dict.subDict("pressureSeries"));
    }
    else if (dict.found("pressureField"))
    {
        Info<< "    pressure is specified as a field" << endl;
        pressureFieldPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    word(dict.lookup("pressureField")),
                    patch().boundaryMesh().mesh().time().timeName(),
                    patch().boundaryMesh().mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh()
            )
        );
    }
    else
    {
        pressure_ = scalarField("pressure", dict, p.size());
    }

    if (secondOrder_)
    {
        Info<< "    second order correction" << endl;
    }

    if (setEffectiveTraction_)
    {
        Info<< "    set effective traction" << endl;
    }

    if (relaxFac_ < 1.0)
    {
        Info<< "    relaxation factor: " << relaxFac_ << endl;
    }
}


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const solidTractionFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(pvf, p, iF, mapper),
    nonOrthogonalCorrections_(pvf.nonOrthogonalCorrections_),
#ifdef OPENFOAM_ORG
    traction_(mapper(pvf.traction_)),
    pressure_(mapper(pvf.pressure_)),
#else
    traction_(pvf.traction_, mapper),
    pressure_(pvf.pressure_, mapper),
#endif
    tractionSeries_(pvf.tractionSeries_),
    pressureSeries_(pvf.pressureSeries_),
    tractionFieldPtr_(),
    pressureFieldPtr_(),
    secondOrder_(pvf.secondOrder_),
    setEffectiveTraction_(pvf.setEffectiveTraction_),
    relaxFac_(pvf.relaxFac_),
    curTimeIndex_(pvf.curTimeIndex_)
{}

#ifndef OPENFOAM_ORG
solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const solidTractionFvPatchVectorField& pvf
)
:
    fixedGradientFvPatchVectorField(pvf),
    nonOrthogonalCorrections_(pvf.nonOrthogonalCorrections_),
    traction_(pvf.traction_),
    pressure_(pvf.pressure_),
    tractionSeries_(pvf.tractionSeries_),
    pressureSeries_(pvf.pressureSeries_),
    tractionFieldPtr_(),
    pressureFieldPtr_(),
    secondOrder_(pvf.secondOrder_),
    setEffectiveTraction_(pvf.setEffectiveTraction_),
    relaxFac_(pvf.relaxFac_),
    curTimeIndex_(pvf.curTimeIndex_)
{}
#endif

solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const solidTractionFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(pvf, iF),
    nonOrthogonalCorrections_(pvf.nonOrthogonalCorrections_),
    traction_(pvf.traction_),
    pressure_(pvf.pressure_),
    tractionSeries_(pvf.tractionSeries_),
    pressureSeries_(pvf.pressureSeries_),
    tractionFieldPtr_(),
    pressureFieldPtr_(),
    secondOrder_(pvf.secondOrder_),
    setEffectiveTraction_(pvf.setEffectiveTraction_),
    relaxFac_(pvf.relaxFac_),
    curTimeIndex_(pvf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);

#ifdef OPENFOAM_ORG
    m(traction_, traction_);
    m(pressure_, pressure_);
#else
    traction_.autoMap(m);
    pressure_.autoMap(m);
#endif
}


// Reverse-map the given fvPatchField onto this fvPatchField
void solidTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& pvf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(pvf, addr);

    const solidTractionFvPatchVectorField& rpvf =
        refCast<const solidTractionFvPatchVectorField>(pvf);

    traction_.rmap(rpvf.traction_, addr);
    pressure_.rmap(rpvf.pressure_, addr);
}


// Update the coefficients associated with the patch field
void solidTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (curTimeIndex_ != db().time().timeIndex())
    {
        curTimeIndex_ = db().time().timeIndex();

        // Called once per time-step

        if (pressureFieldPtr_.valid())
        {
            // Force the pressure field boundary conditions to update
            const_cast<volScalarField&>
            (
                pressureFieldPtr_()
            ).correctBoundaryConditions();
        }

        if (tractionFieldPtr_.valid())
        {
            // Force the traction field boundary conditions to update
            const_cast<volVectorField&>
            (
                tractionFieldPtr_()
            ).correctBoundaryConditions();
        }
    }

    if (tractionFieldPtr_.valid())
    {
        traction_ = tractionFieldPtr_().boundaryField()[patch().index()];
    }
    else if (tractionSeries_.size())
    {
        traction_ = tractionSeries_(this->db().time().timeOutputValue());
    }

    if (pressureFieldPtr_.valid())
    {
        pressure_ = pressureFieldPtr_().boundaryField()[patch().index()];
    }
    else if (pressureSeries_.size())
    {
        pressure_ = pressureSeries_(this->db().time().timeOutputValue());
    }

    scalarField press(pressure_);
    if (setEffectiveTraction_)
    {
        const fvPatchField<scalar>& p =
            patch().lookupPatchField<volScalarField, scalar>("p");

        // Remove the dynamic pressure component: this will force the effective
        // traction to be enforced rather than the total traction
        press -= p;
    }

    // Lookup the solidModel object
    const solidModel& solMod = lookupSolidModel(patch().boundaryMesh().mesh());

    // Set surface-normal gradient on the patch corresponding to the desired
    // traction
    gradient() =
        relaxFac_*solMod.tractionBoundarySnGrad
        (
            traction_, press, patch()
        )
      + (1.0 - relaxFac_)*gradient();

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void solidTractionFvPatchVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    if (nonOrthogonalCorrections_)
    {
        // Lookup the gradient field
        const fvPatchField<tensor>& gradField =
            patch().lookupPatchField<volTensorField, tensor>
            (
            #ifdef OPENFOAM_NOT_EXTEND
                "grad(" + internalField().name() + ")"
            #else
                "grad(" + dimensionedInternalField().name() + ")"
            #endif
            );

        // Non-orthogonal correction vectors
        const vectorField k(patchCorrectionVectors(patch()));

        if (secondOrder_)
        {
            // Face unit normals
            const vectorField n(patch().nf());

            // Correction to internal field cells
            const vectorField dUP(k & gradField.patchInternalField());

            // Normal gradient at internal field cells
            const vectorField nGradUP(n & gradField.patchInternalField());

            Field<vector>::operator=
            (
                patchInternalField()
              + dUP
              + 0.5*(gradient() + nGradUP)/patch().deltaCoeffs()
            );
        }
        else
        {

            Field<vector>::operator=
            (
                patchInternalField()
              + (k & gradField.patchInternalField())
              + gradient()/patch().deltaCoeffs()
            );
        }
    }
    else
    {
        // No non-orthogonal correction
        Field<vector>::operator=
        (
            patchInternalField() + gradient()/patch().deltaCoeffs()
        );
    }

    fvPatchField<vector>::evaluate();
}


void solidTractionFvPatchVectorField::write(Ostream& os) const
{
    // Bug-fix: courtesy of Michael@UW at https://www.cfd-online.com/Forums/
    // openfoam-cc-toolkits-fluid-structure-interaction/221892-solved-paraview
    // -cant-read-solids-files-duplicate-entries-keyword-value.html#post762325
    //fixedGradientFvPatchVectorField::write(os);
    fvPatchVectorField::write(os);

    os.writeKeyword("nonOrthogonalCorrections")
        << nonOrthogonalCorrections_ << token::END_STATEMENT << nl;

    if (tractionFieldPtr_.valid())
    {
        os.writeKeyword("tractionField")
            << tractionFieldPtr_().name() << token::END_STATEMENT << nl;
    }
    else if (tractionSeries_.size())
    {
        os.writeKeyword("tractionSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        tractionSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
#ifdef OPENFOAM_ORG
        writeEntry(os, "traction", traction_);
#else
        traction_.writeEntry("traction", os);
#endif
    }

    if (pressureFieldPtr_.valid())
    {
        os.writeKeyword("pressureField")
            << pressureFieldPtr_().name() << token::END_STATEMENT << nl;
    }
    else if (pressureSeries_.size())
    {
        os.writeKeyword("pressureSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        pressureSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
#ifdef OPENFOAM_ORG
        writeEntry(os, "pressure", pressure_);
#else
        pressure_.writeEntry("pressure", os);
#endif
    }

    os.writeKeyword("secondOrder")
        << secondOrder_ << token::END_STATEMENT << nl;
    os.writeKeyword("setEffectiveTraction")
        << setEffectiveTraction_ << token::END_STATEMENT << nl;
    os.writeKeyword("relaxationFactor")
        << relaxFac_ << token::END_STATEMENT << nl;

#ifdef OPENFOAM_ORG
    writeEntry(os, "value", *this);
    writeEntry(os, "gradient", gradient());
#else
    writeEntry("value", os);
    gradient().writeEntry("gradient", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidTractionFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

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

Class
    Foam::paraboloidInletVelocityFvPatchVectorField

Description
    Sets a paraboloid inlet velocity profile, where there is the option for the
    profile to vary in time following a sine wave.

    Currently only implemented for the x direction.

    Based on the setInletVelocity functionObject by Zeljko Tukovic.

SourceFiles
    paraboloidInletVelocityFvPatchVectorField.C

Author
    Philip Cardiff, UCD. All rights reserved
    Zeljko Tukovic, FSB Zagreb. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef paraboloidInletVelocityFvPatchVectorField_H
#define paraboloidInletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class paraboloidInletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class paraboloidInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Maximum velocity
        const scalar Umax_;

        //- Max y coordinate of paraboloid
        const scalar yMax_;

        //- Max z coordinate of paraboloid
        const scalar zMax_;

        //- Time at which time-varying velocity ends
        //  Set to zero for steady-state velocity
        //  Set to timeAtMaxVelocity for smooth ramp to max velocity
        const scalar timeVaryingEndTime_;

        //- Time that the maximum velocity is reached
        const scalar timeAtMaxVelocity_;

public:

    //- Runtime type information
    TypeName("paraboloidInletVelocity");


    // Constructors

        //- Construct from patch and internal field
        paraboloidInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        paraboloidInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given paraboloidInletVelocityFvPatchVectorField
        //  onto a new patch
        paraboloidInletVelocityFvPatchVectorField
        (
            const paraboloidInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        paraboloidInletVelocityFvPatchVectorField
        (
            const paraboloidInletVelocityFvPatchVectorField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new paraboloidInletVelocityFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        paraboloidInletVelocityFvPatchVectorField
        (
            const paraboloidInletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new paraboloidInletVelocityFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

    // Evaluation functions

        //- Return patch-normal gradient
        virtual tmp<Field<vector> > snGrad() const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<vector> > gradientBoundaryCoeffs() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

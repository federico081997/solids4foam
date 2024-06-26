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
    thermalModel

Description
    Thermal material properties for solids.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    thermalModel.C

\*---------------------------------------------------------------------------*/

#ifndef thermalModel_H
#define thermalModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "thermalLaw.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermalModel Declaration
\*---------------------------------------------------------------------------*/

class thermalModel
:
    public IOdictionary
{
    // Private data

        //- Const reference to the mesh
        const fvMesh& mesh_;

        //- Thermal law
        autoPtr<thermalLaw> lawPtr_;


    // Private Member Functions

        //- Disallow copy construct
        thermalModel(const thermalModel&);

        //- Disallow default bitwise assignment
        void operator=(const thermalModel&);


public:

    //- Runtime type information
    TypeName("thermalModel");

    // Constructors

        //- Construct from dictionary
        thermalModel(const fvMesh& mesh);


    // Destructor

        virtual ~thermalModel()
        {}


    // Member Functions

        //- Return const reference to the mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return thermal law
        const thermalLaw& law() const
        {
            return lawPtr_();
        }

        //- Return specific heat capacity
        tmp<volScalarField> C() const
        {
            return lawPtr_->C();
        }

        //- Return thermal conductivity
        tmp<volScalarField> k() const
        {
            return lawPtr_->k();
        }

        //- Read thermalProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

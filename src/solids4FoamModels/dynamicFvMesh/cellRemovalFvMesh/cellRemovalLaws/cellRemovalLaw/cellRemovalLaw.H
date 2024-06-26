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
    cellRemovalLaw

Description
    Law that selects which cells should be removed

Authors
    Philip Cardif, UCD. All rights reserved.

SourceFiles
    cellRemovalLaw.C
    newCellRemovalLaw.C

\*---------------------------------------------------------------------------*/

#ifndef cellRemovalLaw_H
#define cellRemovalLaw_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class cellRemovalLaw Declaration
\*---------------------------------------------------------------------------*/

class cellRemovalLaw
{
    // Private data

        //- Name
        const word name_;

        //- Reference to stress field
        fvMesh& mesh_;


    // Private Member Functions

        //- Disallow copy construct
        cellRemovalLaw(const cellRemovalLaw&);

        //- Disallow default bitwise assignment
        void operator=(const cellRemovalLaw&);


protected:

        //- Return reference to mesh
        fvMesh& mesh()
        {
            return mesh_;
        }

        //- Return const reference to mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

public:

    //- Runtime type information
    TypeName("cellRemovalLaw");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            cellRemovalLaw,
            dictionary,
            (
                const word name,
                fvMesh& mesh,
                const dictionary& dict
            ),
            (name, mesh, dict)
        );


    // Selectors

        //- Return a reference to the selected meshFailure model
        static autoPtr<cellRemovalLaw> New
        (
            const word& name,
            fvMesh& mesh,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        cellRemovalLaw
        (
            const word& name,
            fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~cellRemovalLaw()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        // Return list of cell IDs to be removed
        virtual tmp<labelField> cellsToRemove() = 0;

        // Return list of patch indices for the faces to be removed
        virtual label exposedFacesPatchID() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

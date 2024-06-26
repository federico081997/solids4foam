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
    solidPointDisplacement

Description
    FunctionObject reports displacement vector value at closest mesh point to
    the specified point.

    The temperature volField is interpolated to a pointField using the
    least squares vol-to-point interpolator.

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    solidPointDisplacement.C

\*---------------------------------------------------------------------------*/

#ifndef solidPointDisplacement_H
#define solidPointDisplacement_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class solidPointDisplacement Declaration
\*---------------------------------------------------------------------------*/

class solidPointDisplacement
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- Optional: region name
        //  This can be "solid" for FSI cases or it can be the name of the
        //  material for multi-material cases
        //  If it is not specified, the "solid" region will still be looked for
        //  to preserve backwards compatibility
        word region_;

        //- Point ID
        label pointID_;

        //- History file ptr
        autoPtr<OFstream> historyFilePtr_;

    // Private Member Functions

        //- Write data
        bool writeData();

        //- Disallow default bitwise copy construct
        solidPointDisplacement
        (
            const solidPointDisplacement&
        );

        //- Disallow default bitwise assignment
        void operator=(const solidPointDisplacement&);


public:

    //- Runtime type information
    TypeName("solidPointDisplacement");


    // Constructors

        //- Construct from components
        solidPointDisplacement
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
#if FOAMEXTEND
        virtual bool execute(const bool forceWrite);
#else
        virtual bool execute();
#endif

        //- Called when time was set at the end of the Time::operator++
        virtual bool timeSet()
        {
            return true;
        }

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);

#ifdef OPENFOAM_NOT_EXTEND
        //- Write
        virtual bool write();
#endif

#ifdef FOAMEXTEND
        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
#endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

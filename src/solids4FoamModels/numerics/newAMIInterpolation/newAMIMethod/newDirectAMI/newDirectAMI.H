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
    Foam::newDirectAMI

Description
    Direct mapped Arbitrary Mesh Interface (AMI) method

SourceFiles
    newDirectAMI.C

\*---------------------------------------------------------------------------*/

#ifdef OPENFOAM_ORG

#ifndef newDirectAMI_H
#define newDirectAMI_H

#include "newAMIMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class newDirectAMI Declaration
\*---------------------------------------------------------------------------*/

template<class SourcePatch, class TargetPatch>
class newDirectAMI
:
    public newAMIMethod<SourcePatch, TargetPatch>
{
    // Private Member Functions

        // Marching front

            //- Append to list of src face seed indices
            void appendToDirectSeeds
            (
                labelList& mapFlag,
                labelList& srcTgtSeed,
                DynamicList<label>& srcSeeds,
                DynamicList<label>& nonOverlapFaces,
                label& srcFacei,
                label& tgtFacei
            ) const;

            //- Restart the advancing front - typically happens for
            //  disconnected regions
            void restartAdvancingFront
            (
                labelList& mapFlag,
                DynamicList<label>& nonOverlapFaces,
                label& srcFacei,
                label& tgtFacei
            ) const;


        // Evaluation

            //- Area of intersection between source and target faces
            scalar interArea
            (
                const label srcFacei,
                const label tgtFacei
            ) const;


public:

    //- Runtime type information
    TypeName("newDirectAMI");


    // Constructors

        //- Construct from components
        newDirectAMI
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const scalarField& srcMagSf,
            const scalarField& tgtMagSf,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool reverseTarget = false,
            const bool requireMatch = true
        );

        //- Disallow default bitwise copy construction
        newDirectAMI(const newDirectAMI&) = delete;


    //- Destructor
    virtual ~newDirectAMI();


    // Member Functions

        // Manipulation

            //- Update addressing and weights
            virtual void calculate
            (
                labelListList& srcAddress,
                scalarListList& srcWeights,
                labelListList& tgtAddress,
                scalarListList& tgtWeights,
                label srcFacei = -1,
                label tgtFacei = -1
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const newDirectAMI&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "newDirectAMI.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

#endif // end of #ifdef OPENFOAM_ORG

// ************************************************************************* //

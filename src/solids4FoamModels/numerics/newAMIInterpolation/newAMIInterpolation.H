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
    Foam::newAMIInterpolation

Description
    Interpolation class dealing with transfer of data between two
    primitive patches with an arbitrary mesh interface (AMI).

    Based on the algorithm given in:

        Conservative interpolation between volume meshes by local Galerkin
        projection, Farrell PE and Maddison JR, 2011, Comput. Methods Appl.
        Mech Engrg, Volume 200, Issues 1-4, pp 89-100

    Interpolation requires that the two patches should have opposite
    orientations (opposite normals).  The 'reverseTarget' flag can be used to
    reverse the orientation of the target patch.


SourceFiles
    newAMIInterpolation.C
    newAMIInterpolationName.C
    newAMIInterpolationParallelOps.C

\*---------------------------------------------------------------------------*/

#ifdef OPENFOAM_ORG

#ifndef newAMIInterpolation_H
#define newAMIInterpolation_H

#include "className.H"
#include "searchableSurface.H"
#include "treeBoundBoxList.H"
#include "boolList.H"
#include "primitivePatch.H"
#include "faceAreaIntersect.H"
#include "globalIndex.H"
#include "ops.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class newAMIInterpolationName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(newAMIInterpolation);


/*---------------------------------------------------------------------------*\
                      Class newAMIInterpolation Declaration
\*---------------------------------------------------------------------------*/

template<class SourcePatch, class TargetPatch>
class newAMIInterpolation
:
    public newAMIInterpolationName
{
public:

    // Public typedefs

        // Typedef to SourcePatch type this newAMIInterplation is instantiated on
        typedef SourcePatch sourcePatchType;

        // Typedef to TargetPatch type this newAMIInterplation is instantiated on
        typedef TargetPatch targetPatchType;


    // Public data types

        //- Enumeration specifying interpolation method
        enum interpolationMethod
        {
            imDirect,
            imMapNearest,
            imFaceAreaWeight,
            imPartialFaceAreaWeight,
            imSweptFaceAreaWeight
        };

        //- Convert interpolationMethod to word representation
        static word interpolationMethodToWord
        (
            const interpolationMethod& method
        );

        //- Convert word to interpolationMethod
        static interpolationMethod wordTointerpolationMethod
        (
            const word& method
        );

        //- Calculate the patch face magnitudes for the given tri-mode
        template<class Patch>
        static tmp<scalarField> patchMagSf
        (
            const Patch& patch,
            const faceAreaIntersect::triangulationMode triMode
        );


private:

    // Private Data

        //- Interpolation method
        const word methodName_;

        //- Flag to indicate that the two patches are co-directional and
        //  that the orientation of the target patch should be reversed
        const bool reverseTarget_;

        //- Flag to indicate that the two patches must be matched/an overlap
        //  exists between them
        const bool requireMatch_;

        //- Index of processor that holds all of both sides. -1 in all other
        //  cases
        label singlePatchProc_;

        //- Threshold weight below which interpolation is deactivated
        scalar lowWeightCorrection_;


        // Source patch

            //- Source face areas
            scalarField srcMagSf_;

            //- Addresses of target faces per source face
            labelListList srcAddress_;

            //- Weights of target faces per source face
            scalarListList srcWeights_;

            //- Sum of weights of target faces per source face
            scalarField srcWeightsSum_;


        // Target patch

            //- Target face areas
            scalarField tgtMagSf_;

            //- Addresses of source faces per target face
            labelListList tgtAddress_;

            //- Weights of source faces per target face
            scalarListList tgtWeights_;

            //- Sum of weights of source faces per target face
            scalarField tgtWeightsSum_;


        //- Face triangulation mode
        const faceAreaIntersect::triangulationMode triMode_;

        //- Source map pointer - parallel running only
        autoPtr<mapDistribute> srcMapPtr_;

        //- Target map pointer - parallel running only
        autoPtr<mapDistribute> tgtMapPtr_;

        //- Patches are globalPolyPatches
        //  same copy of global patches on all processors
        const bool useGlobalPolyPatch_;


        // Parallel functionality

            //- Calculate if patches are on multiple processors
            label calcDistribution
            (
                const SourcePatch& srcPatch,
                const TargetPatch& tgtPatch
            ) const;

            label calcOverlappingProcs
            (
                const List<treeBoundBoxList>& procBb,
                const treeBoundBox& bb,
                boolList& overlaps
            ) const;

            void distributePatches
            (
                const mapDistribute& map,
                const TargetPatch& pp,
                const globalIndex& gi,
                List<faceList>& faces,
                List<pointField>& points,
                List<labelList>& tgtFaceIDs
            ) const;

            void distributeAndMergePatches
            (
                const mapDistribute& map,
                const TargetPatch& tgtPatch,
                const globalIndex& gi,
                faceList& tgtFaces,
                pointField& tgtPoints,
                labelList& tgtFaceIDs
            ) const;

            autoPtr<mapDistribute> calcProcMap
            (
                const SourcePatch& srcPatch,
                const TargetPatch& tgtPatch
            ) const;


        // Initialisation

            //- Project points to surface
            void projectPointsToSurface
            (
                const searchableSurface& surf,
                pointField& pts
            ) const;


        // Manipulation

            //- Sum the weights for each face
            static void sumWeights
            (
                const scalarListList& wght,
                scalarField& wghtSum
            );

            //- As above, but for multiple sets of weights
            static void sumWeights
            (
                const UPtrList<scalarListList>& wghts,
                scalarField& wghtSum
            );

            //- Print out information relating to the weights sum. Values close
            //  to one are ideal. This information acts as a measure of the
            //  quality of the AMI.
            static void reportSumWeights
            (
                const scalarField& patchAreas,
                const word& patchName,
                const scalarField& wghtSum,
                const scalar lowWeightTol
            );

            //- Normalise the weights so that they sum to one for each face.
            //  This may stabilise the solution at the expense of accuracy.
            static void normaliseWeights
            (
                scalarListList& wght,
                const scalarField& wghtSum
            );

            //- As above but for multiple sets of weights
            static void normaliseWeights
            (
                UPtrList<scalarListList>& wghts,
                const scalarField& wghtSum
            );


        // Constructor helpers

            static void agglomerate
            (
                const autoPtr<mapDistribute>& targetMap,
                const scalarField& fineSrcMagSf,
                const labelListList& fineSrcAddress,
                const scalarListList& fineSrcWeights,

                const labelList& sourceRestrictAddressing,
                const labelList& targetRestrictAddressing,

                scalarField& srcMagSf,
                labelListList& srcAddress,
                scalarListList& srcWeights,
                scalarField& srcWeightsSum,
                autoPtr<mapDistribute>& tgtMap
            );

            void constructFromSurface
            (
                const SourcePatch& srcPatch,
                const TargetPatch& tgtPatch,
                const autoPtr<searchableSurface>& surfPtr,
                const bool report
            );

public:

    // Constructors

        //- Construct from components
        newAMIInterpolation
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool requireMatch = true,
            const interpolationMethod& method = imFaceAreaWeight,
            const scalar lowWeightCorrection = -1,
            const bool reverseTarget = false,
            const bool useGlobalPolyPatch = false,
            const bool report = true
        );

        //- Construct from components
        newAMIInterpolation
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool requireMatch = true,
            const word& methodName =
                interpolationMethodToWord(imFaceAreaWeight),
            const scalar lowWeightCorrection = -1,
            const bool reverseTarget = false,
            const bool useGlobalPolyPatch = false,
            const bool report = true
        );

        //- Construct from components, with projection surface
        newAMIInterpolation
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const autoPtr<searchableSurface>& surf,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool requireMatch = true,
            const interpolationMethod& method = imFaceAreaWeight,
            const scalar lowWeightCorrection = -1,
            const bool reverseTarget = false,
            const bool useGlobalPolyPatch = false,
            const bool report = true
        );

        //- Construct from components, with projection surface
        newAMIInterpolation
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const autoPtr<searchableSurface>& surf,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool requireMatch = true,
            const word& methodName =
                interpolationMethodToWord(imFaceAreaWeight),
            const scalar lowWeightCorrection = -1,
            const bool reverseTarget = false,
            const bool useGlobalPolyPatch = false,
            const bool report = true
        );

        //- Construct from agglomeration of newAMIInterpolation.
        //  Agglomeration passed in as new coarse size and addressing from fine
        //  from coarse
        newAMIInterpolation
        (
            const newAMIInterpolation<SourcePatch, TargetPatch>& fineAMI,
            const labelList& sourceRestrictAddressing,
            const labelList& neighbourRestrictAddressing,
            const bool report = false
        );

        //- Disallow default bitwise copy construction
        newAMIInterpolation(const newAMIInterpolation&) = delete;


    //- Destructor
    ~newAMIInterpolation();


    // Member Functions

        // Access

            //- Set to -1, or the processor holding all faces (both sides) of
            //  the AMI
            inline label singlePatchProc() const;

            //- Threshold weight below which interpolation is deactivated
            inline scalar lowWeightCorrection() const;

            //- Return true if employing a 'lowWeightCorrection'
            inline bool applyLowWeightCorrection() const;


            // Source patch

                //- Return const access to source patch face areas
                inline const scalarField& srcMagSf() const;

                //- Return const access to source patch addressing
                inline const labelListList& srcAddress() const;

                //- Return const access to source patch weights
                inline const scalarListList& srcWeights() const;

                //- Return access to source patch weights
                inline scalarListList& srcWeights();

                //- Return const access to normalisation factor of source
                //  patch weights (i.e. the sum before normalisation)
                inline const scalarField& srcWeightsSum() const;

                //- Return access to normalisation factor of source
                //  patch weights (i.e. the sum before normalisation)
                inline scalarField& srcWeightsSum();

                //- Source map pointer - valid only if singlePatchProc = -1
                //  This gets source data into a form to be consumed by
                //  tgtAddress, tgtWeights
                inline const mapDistribute& srcMap() const;


            // Target patch

                //- Return const access to target patch face areas
                inline const scalarField& tgtMagSf() const;

                //- Return const access to target patch addressing
                inline const labelListList& tgtAddress() const;

                //- Return const access to target patch weights
                inline const scalarListList& tgtWeights() const;

                //- Return access to target patch weights
                inline scalarListList& tgtWeights();

                //- Return const access to normalisation factor of target
                //  patch weights (i.e. the sum before normalisation)
                inline const scalarField& tgtWeightsSum() const;

                //- Return access to normalisation factor of target
                //  patch weights (i.e. the sum before normalisation)
                inline scalarField& tgtWeightsSum();

                //- Target map pointer -  valid only if singlePatchProc=-1.
                //  This gets target data into a form to be consumed by
                //  srcAddress, srcWeights
                inline const mapDistribute& tgtMap() const;


        // Manipulation

            //- Update addressing and weights
            void update
            (
                const SourcePatch& srcPatch,
                const TargetPatch& tgtPatch,
                const bool report
            );

            //- Sum the weights on both sides of an AMI
            static void sumWeights
            (
                newAMIInterpolation<SourcePatch, TargetPatch>& AMI
            );

            //- As above, but for multiple AMI-s
            static void sumWeights
            (
                PtrList<newAMIInterpolation<SourcePatch, TargetPatch>>& AMIs
            );

            //- Print out information relating to the weights sum. Values close
            //  to one are ideal. This information acts as a measure of the
            //  quality of the AMI.
            static void reportSumWeights
            (
                newAMIInterpolation<SourcePatch, TargetPatch>& AMI
            );

            //- Normalise the weights on both sides of an AMI
            static void normaliseWeights
            (
                newAMIInterpolation<SourcePatch, TargetPatch>& AMI
            );

            //- As above, but for multiple AMI-s
            static void normaliseWeights
            (
                UPtrList<newAMIInterpolation<SourcePatch, TargetPatch>>& AMIs
            );


        // Evaluation

            // Low-level

                //- Interpolate from target to source with supplied op
                //  to combine existing value with remote value and weight
                template<class Type, class CombineOp>
                void interpolateToSource
                (
                    const UList<Type>& fld,
                    const CombineOp& cop,
                    List<Type>& result,
                    const UList<Type>& defaultValues = UList<Type>::null()
                ) const;

                //- Interpolate from source to target with supplied op
                //  to combine existing value with remote value and weight
                template<class Type, class CombineOp>
                void interpolateToTarget
                (
                    const UList<Type>& fld,
                    const CombineOp& cop,
                    List<Type>& result,
                    const UList<Type>& defaultValues = UList<Type>::null()
                ) const;


            //- Interpolate from target to source with supplied op
            template<class Type, class CombineOp>
            tmp<Field<Type>> interpolateToSource
            (
                const Field<Type>& fld,
                const CombineOp& cop,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from target tmp field to source with supplied op
            template<class Type, class CombineOp>
            tmp<Field<Type>> interpolateToSource
            (
                const tmp<Field<Type>>& tFld,
                const CombineOp& cop,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from source to target with supplied op
            template<class Type, class CombineOp>
            tmp<Field<Type>> interpolateToTarget
            (
                const Field<Type>& fld,
                const CombineOp& cop,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from source tmp field to target with supplied op
            template<class Type, class CombineOp>
            tmp<Field<Type>> interpolateToTarget
            (
                const tmp<Field<Type>>& tFld,
                const CombineOp& cop,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from target to source
            template<class Type>
            tmp<Field<Type>> interpolateToSource
            (
                const Field<Type>& fld,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from target tmp field
            template<class Type>
            tmp<Field<Type>> interpolateToSource
            (
                const tmp<Field<Type>>& tFld,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from source to target
            template<class Type>
            tmp<Field<Type>> interpolateToTarget
            (
                const Field<Type>& fld,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from source tmp field
            template<class Type>
            tmp<Field<Type>> interpolateToTarget
            (
                const tmp<Field<Type>>& tFld,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;


        // Point intersections

            //- Return source patch face index of point on target patch face
            label srcPointFace
            (
                const SourcePatch& srcPatch,
                const TargetPatch& tgtPatch,
                const vector& n,
                const label tgtFacei,
                point& tgtPoint
            )
            const;

            //- Return target patch face index of point on source patch face
            label tgtPointFace
            (
                const SourcePatch& srcPatch,
                const TargetPatch& tgtPatch,
                const vector& n,
                const label srcFacei,
                point& srcPoint
            )
            const;


        // Checks

            //- Write face connectivity as OBJ file
            void writeFaceConnectivity
            (
                const SourcePatch& srcPatch,
                const TargetPatch& tgtPatch,
                const labelListList& srcAddress
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const newAMIInterpolation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "newAMIInterpolationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "newAMIInterpolation.C"
    #include "newAMIInterpolationParallelOps.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

#endif // end of #ifdef OPENFOAM_ORG

// ************************************************************************* //

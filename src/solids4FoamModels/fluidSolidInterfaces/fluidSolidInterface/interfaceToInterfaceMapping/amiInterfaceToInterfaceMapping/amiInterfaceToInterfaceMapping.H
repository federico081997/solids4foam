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
    amiInterfaceToInterfaceMapping

Description
    interfaceToInterfaceMapping wrapper for AMI interpolation

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    amiInterfaceToInterfaceMapping.C

\*---------------------------------------------------------------------------*/

#ifdef OPENFOAM_NOT_EXTEND

#ifndef amiInterfaceToInterfaceMapping_H
#define amiInterfaceToInterfaceMapping_H

#include "interfaceToInterfaceMapping.H"
#include "amiZoneInterpolation.H"
#include "FieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace interfaceToInterfaceMappings
{

/*---------------------------------------------------------------------------*\
                        Class amiInterfaceToInterfaceMapping Declaration
\*---------------------------------------------------------------------------*/

class amiInterfaceToInterfaceMapping
:
    public interfaceToInterfaceMapping
{
    // Private data

            //- AMI interpolator between two standAlone patches
            mutable autoPtr<amiZoneInterpolation> interpolatorPtr_;

            //- zoneB face into which each point of zoneA is projected
            mutable List<labelPair>* zoneAPointAddressingPtr_;

            //- zoneA point weighting factors
            mutable FieldField<Field, scalar>* zoneAPointWeightsPtr_;

            //- zoneA face into which each point of zoneB is projected
            mutable List<labelPair>* zoneBPointAddressingPtr_;

            //- zoneB point weighting factors
            mutable FieldField<Field, scalar>* zoneBPointWeightsPtr_;

    // Private Member Functions

        //- Make the AMI interpolator
        void makeInterpolator() const;

        //- Return a const access to the AMI interpolator
        const amiZoneInterpolation& interpolator() const;

        //- Check zoneA to zoneB face interpolation error
        void checkZoneAToZoneBError() const;

        //- Check zoneB to zoneA point interpolation error
        //void checkZoneBToZoneAError() const;

        //- Calculate zoneA point addressing
        void calcZoneAPointAddressing() const;

        //- Calculate zoneA point weights
        void calcZoneAPointWeights() const;

        //- Calculate zoneB point addressing
        void calcZoneBPointAddressing() const;

        //- Calculate zoneB point weights
        void calcZoneBPointWeights() const;

        //- Return reference to zoneA point addressing
        const List<labelPair>& zoneAPointAddr() const;

        //- Return reference to zoneA point weights
        const FieldField<Field, scalar>& zoneAPointWeights() const;

        //- Return reference to zoneB point addressing
        const List<labelPair>& zoneBPointAddr() const;

        //- Return reference to zoneB point weights
        const FieldField<Field, scalar>& zoneBPointWeights() const;

        //- Transfer/map/interpolate from one zone faces to another zone
        //  faces for Type
        template<class Type>
        void transferFacesZoneToZone
        (
            const standAlonePatch& fromZone,    // from zone
            const standAlonePatch& toZone,      // to zone
            const Field<Type>& fromField,       // from field
            Field<Type>& toField               // to field
        ) const;

        //- Transfer/map/interpolate from one zone points to another zone
        //  points for Type
        template<class Type>
        void transferPointsZoneToZone
        (
            const standAlonePatch& fromZone,    // from zone
            const standAlonePatch& toZone,      // to zone
            const Field<Type>& fromField,       // from field
            Field<Type>& toField               // to field
        ) const;

        //- Disallow default bitwise copy construct
        amiInterfaceToInterfaceMapping(const amiInterfaceToInterfaceMapping&);

        //- Disallow default bitwise assignment
        void operator=(const amiInterfaceToInterfaceMapping&);


public:

    //- Runtime type information
    TypeName("AMI");


    // Constructors

        //- Construct from components
        amiInterfaceToInterfaceMapping
        (
            const word& type,
            const dictionary& dict,
            const primitivePatch& patchA,
            const primitivePatch& patchB,
            const globalPolyPatch& globalPatchA,
            const globalPolyPatch& globalPatchB
        );


    // Destructor

        virtual ~amiInterfaceToInterfaceMapping()
        {}


    // Member Functions

        // Edit

            //- Transfer/map/interpolate from one zone faces to another zone
            //  faces for scalars
            virtual void transferFacesZoneToZone
            (
                const standAlonePatch& fromZone,    // from zone
                const standAlonePatch& toZone,      // to zone
                const Field<scalar>& fromField,     // from field
                Field<scalar>& toField              // to field
            ) const;

            //- Transfer/map/interpolate from one zone points to another zone
            //  points for scalars
            virtual void transferPointsZoneToZone
            (
                const standAlonePatch& fromZone,    // from zone
                const standAlonePatch& toZone,      // to zone
                const Field<scalar>& fromField,     // from field
                Field<scalar>& toField              // to field
            ) const;

            //- Transfer/map/interpolate from one zone faces to another zone
            //  faces for vectors
            virtual void transferFacesZoneToZone
            (
                const standAlonePatch& fromZone,      // from zone
                const standAlonePatch& toZone,        // to zone
                const Field<vector>& fromField,       // from field
                Field<vector>& toField                // to field
            ) const;

            //- Transfer/map/interpolate from one zone points to another zone
            //  points for vectors
            virtual void transferPointsZoneToZone
            (
                const standAlonePatch& fromZone,      // from zone
                const standAlonePatch& toZone,        // to zone
                const Field<vector>& fromField,       // from field
                Field<vector>& toField                // to field
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfaceToInterfaceMappings

} // End namespace Foam

#ifdef NoRepository
#   include "amiInterfaceToInterfaceMappingTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

#endif // end of #ifdef OPENFOAM_NOT_EXTEND

// ************************************************************************* //

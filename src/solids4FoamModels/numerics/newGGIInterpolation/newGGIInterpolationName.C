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

Description
    Shared template name for GGI interpolation

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

Contributor:
    Martin Beaudoin, Hydro-Quebec, (2008)

\*---------------------------------------------------------------------------*/

#include "newGGIInterpolationTemplate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::newGGIInterpolationName, 0);


template<>
const char*
Foam::NamedEnum<Foam::newGGIInterpolationName::quickReject, 3>::names[] =
{
    "distance3D",
    "AABB",
    "bbOctree"
};


const Foam::NamedEnum<Foam::newGGIInterpolationName::quickReject, 3>
    Foam::newGGIInterpolationName::quickRejectNames_;


// ************************************************************************* //

// Created by: Frank J. Mattia
// Copyright (c) 2026 OPEN CASCADE SAS
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.
//
// Alternatively, this file may be used under the terms of Open CASCADE
// commercial license or contractual agreement.

#ifndef _BRepAlgoAPI_SolidSweep_HeaderFile
#define _BRepAlgoAPI_SolidSweep_HeaderFile

#include <NCollection_Sequence.hxx>
#include <Standard.hxx>
#include <TopoDS_Solid.hxx>
#include <gp_Ax1.hxx>
#include <gp_Vec.hxx>

class BRepAlgoAPI_SolidSweep
{
public:
  //! Circular translation around a world-space axis with fixed body orientation.
  struct CircularTranslation
  {
    gp_Ax1 Axis;
    double Angle = 0.0; //!< Signed angle in radians.
  };

  //! Sweeps a valid solid along a straight translation.
  //! Supports planar-faced solids and full right circular cylinders.
  //! Other cases raise Standard_ConstructionError.
  Standard_EXPORT static TopoDS_Solid Sweep(const TopoDS_Solid& theSolid,
                                            const gp_Vec&       theTranslation);

  //! Sweeps a solid along consecutive translation offsets from its original pose.
  //! Each leg supports the same solids as straight translation.
  //! One offset places the solid at that pose. An empty path is invalid.
  Standard_EXPORT static TopoDS_Solid Sweep(const TopoDS_Solid&                 theSolid,
                                            const NCollection_Sequence<gp_Vec>& theOffsets);

  //! Sweeps a solid along a circular arc with fixed orientation.
  //! Supports full right circular cylinders with axis parallel to the path axis
  //! and signed angles in [-pi, pi]. Other cases raise Standard_ConstructionError.
  Standard_EXPORT static TopoDS_Solid Sweep(const TopoDS_Solid&        theSolid,
                                            const CircularTranslation& theMotion);
};

#endif

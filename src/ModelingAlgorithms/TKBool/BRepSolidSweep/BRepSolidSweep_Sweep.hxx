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

#ifndef _BRepSolidSweep_Sweep_HeaderFile
#define _BRepSolidSweep_Sweep_HeaderFile

#include <NCollection_Sequence.hxx>
#include <Standard.hxx>
#include <TopoDS_Solid.hxx>
#include <gp_Ax1.hxx>
#include <gp_Vec.hxx>

class BRepSolidSweep_Sweep
{
public:
  //! Circular translation around a world-space axis with fixed body orientation.
  struct CircularTranslation
  {
    gp_Ax1 Axis;
    double Angle = 0.0; //!< Signed angle in radians.
  };

  //! Helical translation around a world-space axis with fixed body orientation.
  struct HelicalTranslation
  {
    gp_Ax1 Axis;
    double Angle            = 0.0;    //!< Signed angle in radians.
    double AdvancePerTurn   = 0.0;    //!< Axial advance along Axis per revolution of travel.
    double SurfaceTolerance = 1.0e-5; //!< Target sampled deviation in model units.
  };

  //! Rotation about a fixed axis combined with translation along that axis.
  struct ScrewMotion
  {
    gp_Ax1 Axis;
    double Angle            = 0.0;   //!< Positive angle in radians.
    double AdvancePerTurn   = 0.0;   //!< Signed axial advance per revolution.
    double SurfaceTolerance = 0.005; //!< Target sampled deviation in model units.
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
  //! Builds smooth envelope patches for a full right circular cylinder on a helix.
  //! Raises Standard_ConstructionError if a valid swept solid cannot be built.
  Standard_EXPORT static TopoDS_Solid Sweep(const TopoDS_Solid&       theSolid,
                                            const HelicalTranslation& theMotion);
  //! Sweeps an axisymmetric solid with two disk caps and tangent side faces.
  //! Supports coaxial cylindrical, conical and toroidal bands, positive
  //! cap-normal velocity and non-overlapping turns.
  Standard_EXPORT static TopoDS_Solid Sweep(const TopoDS_Solid& theSolid,
                                            const ScrewMotion& theMotion);
};

#endif

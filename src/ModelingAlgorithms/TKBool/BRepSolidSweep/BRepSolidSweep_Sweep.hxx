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
#include <gp_Trsf.hxx>

#include <functional>

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

  //! Rotation about Axis after translation by station * Drift.
  //! The angular rate is InitialRate / (1 + station * RateDenominatorSlope).
  struct RevolvingMotion
  {
    gp_Ax1 Axis;
    gp_Vec Drift;
    double InitialRate          = 0.0;
    double RateDenominatorSlope = 0.0;
    double Length               = 0.0;   //!< Final value of the station parameter.
    double SurfaceTolerance     = 0.005; //!< Target sampled deviation in model units.
  };

  //! Instantaneous rigid placement and its derivatives in world coordinates.
  struct RigidPose
  {
    gp_Trsf Placement;
    gp_Vec  AngularVelocity; //!< World angular velocity per station.
    gp_Vec  OriginVelocity;  //!< Derivative of Placement.TranslationPart().
  };

  //! Smooth motion over a finite station interval. Placement maps the source
  //! solid into world coordinates, including at the initial station.
  struct RigidMotion
  {
    std::function<RigidPose(double)> Evaluate;
    double                           First            = 0.0;
    double                           Last             = 1.0;
    double                           SurfaceTolerance = 0.005;
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
                                            const ScrewMotion&  theMotion);
  //! Sweeps analytic faces and sharp edges whose contact loops persist
  //! continuously. Intersects the candidate envelope and keeps occupied cells.
  //! Raises Standard_ConstructionError if the contact topology changes or the
  //! occupied region cannot be recovered as a single valid solid.
  Standard_EXPORT static TopoDS_Solid Sweep(const TopoDS_Solid&    theSolid,
                                            const RevolvingMotion& theMotion);
  //! Builds and trims a moving contact envelope for an analytic source solid
  //! with persistent contact loops. The supplied derivatives must agree with
  //! the placement law. Raises Standard_ConstructionError on unsupported
  //! contact topology or if the result cannot be closed.
  Standard_EXPORT static TopoDS_Solid Sweep(const TopoDS_Solid& theSolid,
                                             const RigidMotion& theMotion);
};

#endif

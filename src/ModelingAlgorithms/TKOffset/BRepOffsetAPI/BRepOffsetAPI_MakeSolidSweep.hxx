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

#ifndef _BRepOffsetAPI_MakeSolidSweep_HeaderFile
#define _BRepOffsetAPI_MakeSolidSweep_HeaderFile

#include <BRepBuilderAPI_MakeShape.hxx>
#include <BRepSolidSweep_Sweep.hxx>

#include <variant>

//! Builds the volume occupied by a solid undergoing rigid motion.
//! Subshape history is not provided by this operation.
class BRepOffsetAPI_MakeSolidSweep : public BRepBuilderAPI_MakeShape
{
public:
  using CircularTranslation = BRepSolidSweep_Sweep::CircularTranslation;
  using HelicalTranslation  = BRepSolidSweep_Sweep::HelicalTranslation;
  using ScrewMotion         = BRepSolidSweep_Sweep::ScrewMotion;
  using RigidPose           = BRepSolidSweep_Sweep::RigidPose;
  using RigidMotion         = BRepSolidSweep_Sweep::RigidMotion;
  using RevolvingMotion     = BRepSolidSweep_Sweep::RevolvingMotion;
  using Motion              = std::variant<gp_Vec,
                                           NCollection_Sequence<gp_Vec>,
                                           CircularTranslation,
                                           HelicalTranslation,
                                           ScrewMotion>;

  //! Initializes the builder. Build() or Solid() computes the swept volume.
  Standard_EXPORT BRepOffsetAPI_MakeSolidSweep(const TopoDS_Solid& theSolid,
                                               const Motion&       theMotion);

  Standard_EXPORT void Build(
    const Message_ProgressRange& theRange = Message_ProgressRange()) override;

  //! Returns the swept solid, building it if necessary.
  Standard_EXPORT TopoDS_Solid Solid();

  //! Convenience entry point for constructing and returning a swept solid.
  Standard_EXPORT static TopoDS_Solid Sweep(const TopoDS_Solid& theSolid, const Motion& theMotion);

  Standard_EXPORT static TopoDS_Solid Sweep(const TopoDS_Solid&    theSolid,
                                            const RevolvingMotion& theMotion);

  Standard_EXPORT static TopoDS_Solid Sweep(const TopoDS_Solid& theSolid,
                                            const RigidMotion&  theMotion);

private:
  TopoDS_Solid mySolid;
  Motion       myMotion;
};

#endif

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

#include <BRepOffsetAPI_MakeSolidSweep.hxx>

#include <TopoDS.hxx>

BRepOffsetAPI_MakeSolidSweep::BRepOffsetAPI_MakeSolidSweep(const TopoDS_Solid& theSolid,
                                                           const Motion&       theMotion)
    : mySolid(theSolid),
      myMotion(theMotion)
{
}

void BRepOffsetAPI_MakeSolidSweep::Build(const Message_ProgressRange& theRange)
{
  NotDone();
  myShape.Nullify();
  if (theRange.UserBreak())
  {
    return;
  }
  const TopoDS_Solid aResult = std::visit(
    [this](const auto& theMotion) { return BRepSolidSweep_Sweep::Sweep(mySolid, theMotion); },
    myMotion);
  if (theRange.More())
  {
    myShape = aResult;
    Done();
  }
}

TopoDS_Solid BRepOffsetAPI_MakeSolidSweep::Solid()
{
  return TopoDS::Solid(Shape());
}

TopoDS_Solid BRepOffsetAPI_MakeSolidSweep::Sweep(const TopoDS_Solid& theSolid,
                                                 const Motion&       theMotion)
{
  BRepOffsetAPI_MakeSolidSweep aBuilder(theSolid, theMotion);
  return aBuilder.Solid();
}

TopoDS_Solid BRepOffsetAPI_MakeSolidSweep::Sweep(const TopoDS_Solid&    theSolid,
                                                 const RevolvingMotion& theMotion)
{
  return BRepSolidSweep_Sweep::Sweep(theSolid, theMotion);
}

TopoDS_Solid BRepOffsetAPI_MakeSolidSweep::Sweep(const TopoDS_Solid& theSolid,
                                                 const RigidMotion&  theMotion)
{
  return BRepSolidSweep_Sweep::Sweep(theSolid, theMotion);
}

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

#include <BRepAlgoAPI_SolidSweep.hxx>

#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <NCollection_List.hxx>
#include <Precision.hxx>
#include <Standard_ConstructionError.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Pln.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>

#include <algorithm>
#include <cmath>

namespace
{
bool IsFinite(const gp_Vec& theVector)
{
  return std::isfinite(theVector.X()) && std::isfinite(theVector.Y())
         && std::isfinite(theVector.Z());
}

TopoDS_Solid FuseSweep(const TopoDS_Solid&                   theStart,
                       const NCollection_List<TopoDS_Shape>& theTools,
                       bool                                  theSimplify = false)
{
  if (theTools.IsEmpty())
  {
    return theStart;
  }

  NCollection_List<TopoDS_Shape> anArguments;
  anArguments.Append(theStart);
  BRepAlgoAPI_Fuse aFuse;
  aFuse.SetArguments(anArguments);
  aFuse.SetTools(theTools);
  aFuse.SetNonDestructive(theSimplify);
  aFuse.Build();
  if (aFuse.HasErrors())
  {
    throw Standard_ConstructionError("Unable to fuse the swept volumes");
  }
  if (theSimplify)
  {
    aFuse.SimplifyResult();
  }

  TopExp_Explorer aSolids(aFuse.Shape(), TopAbs_SOLID);
  if (!aSolids.More())
  {
    throw Standard_ConstructionError("The sweep produced no solid");
  }
  const TopoDS_Solid aResult = TopoDS::Solid(aSolids.Current());
  aSolids.Next();
  if (aSolids.More())
  {
    throw Standard_ConstructionError("The sweep produced more than one solid");
  }
  return aResult;
}

TopoDS_Solid SolidAtOffset(const TopoDS_Solid& theSolid, const gp_Vec& theOffset)
{
  gp_Trsf aPlacement;
  aPlacement.SetTranslation(theOffset);
  return TopoDS::Solid(theSolid.Moved(TopLoc_Location(aPlacement)));
}

struct CylinderGeometry
{
  gp_Pnt      Base;
  gp_Dir      Axis;
  double      Radius;
  double      Height;
  TopoDS_Face LowerCap;
  TopoDS_Face UpperCap;
};

CylinderGeometry ReadCylinder(const TopoDS_Solid& theCylinder);
TopoDS_Solid     SweepCylinderTransversely(const TopoDS_Solid& theCylinder,
                                           const gp_Vec&       theTranslation);
} // namespace

TopoDS_Solid BRepAlgoAPI_SolidSweep::Sweep(const TopoDS_Solid& theSolid,
                                           const gp_Vec&       theTranslation)
{
  if (theSolid.IsNull() || !IsFinite(theTranslation))
  {
    throw Standard_ConstructionError("Expected a solid and a finite translation");
  }

  NCollection_List<TopoDS_Shape> aTools;

  for (TopExp_Explorer anExplorer(theSolid, TopAbs_FACE); anExplorer.More(); anExplorer.Next())
  {
    const TopoDS_Face&        aFace = TopoDS::Face(anExplorer.Current());
    const BRepAdaptor_Surface aSurface(aFace);
    if (aSurface.GetType() == GeomAbs_Cylinder)
    {
      const gp_Vec anAxis(aSurface.Cylinder().Axis().Direction());
      const double anAxialDistance = theTranslation.Dot(anAxis);
      if (std::abs(anAxialDistance) <= Precision::Confusion())
      {
        return SweepCylinderTransversely(theSolid, theTranslation);
      }
      const gp_Vec           anAxialMotion     = anAxialDistance * anAxis;
      const gp_Vec           aTransverseMotion = theTranslation - anAxialMotion;
      const CylinderGeometry aCylinder         = ReadCylinder(theSolid);
      if (aTransverseMotion.Magnitude() > Precision::Confusion())
      {
        // In this frame, the outward normals on the leading half of the
        // cylinder have positive dot product with the transverse motion.
        const gp_Ax3            aFrame(aCylinder.Base, aCylinder.Axis, gp_Dir(aTransverseMotion));
        const double            aHalfPi = 0.5 * std::acos(-1.0);
        BRepBuilderAPI_MakeFace aLateralFace(gp_Cylinder(aFrame, aCylinder.Radius),
                                             -aHalfPi,
                                             aHalfPi,
                                             0.0,
                                             aCylinder.Height);
        const TopoDS_Face& aCap = anAxialDistance > 0.0 ? aCylinder.UpperCap : aCylinder.LowerCap;
        aTools.Clear();
        aTools.Append(BRepPrimAPI_MakePrism(aLateralFace.Face(), theTranslation, true).Shape());
        aTools.Append(BRepPrimAPI_MakePrism(aCap, theTranslation, true).Shape());
        return FuseSweep(theSolid, aTools);
      }

      const gp_Pnt aBase = aCylinder.Base.Translated(std::min(0.0, anAxialDistance) * anAxis);
      return BRepPrimAPI_MakeCylinder(gp_Ax2(aBase, aCylinder.Axis),
                                      aCylinder.Radius,
                                      aCylinder.Height + std::abs(anAxialDistance))
        .Solid();
    }
    if (aSurface.GetType() != GeomAbs_Plane)
    {
      throw Standard_ConstructionError("Unsupported face type for solid translation");
    }

    gp_Dir aNormal = aSurface.Plane().Axis().Direction();
    if (aFace.Orientation() == TopAbs_REVERSED)
    {
      aNormal.Reverse();
    }

    if (gp_Vec(aNormal).Dot(theTranslation) > Precision::Confusion())
    {
      BRepPrimAPI_MakePrism aPrism(aFace, theTranslation, true);
      aTools.Append(aPrism.Shape());
    }
  }

  return FuseSweep(theSolid, aTools);
}

TopoDS_Solid BRepAlgoAPI_SolidSweep::Sweep(const TopoDS_Solid&                 theSolid,
                                           const NCollection_Sequence<gp_Vec>& theOffsets)
{
  if (theSolid.IsNull() || theOffsets.IsEmpty())
  {
    throw Standard_ConstructionError("Expected a solid and at least one offset");
  }
  for (const gp_Vec& anOffset : theOffsets)
  {
    if (!IsFinite(anOffset))
    {
      throw Standard_ConstructionError("Expected finite polyline offsets");
    }
  }

  const TopoDS_Solid             aStart = SolidAtOffset(theSolid, theOffsets.First());
  NCollection_List<TopoDS_Shape> aLegs;
  for (int anIndex = 2; anIndex <= theOffsets.Length(); ++anIndex)
  {
    const gp_Vec& aFrom   = theOffsets.Value(anIndex - 1);
    const gp_Vec  aMotion = theOffsets.Value(anIndex) - aFrom;
    if (aMotion.SquareMagnitude() == 0.0)
    {
      continue;
    }
    aLegs.Append(Sweep(SolidAtOffset(theSolid, aFrom), aMotion));
  }
  // Merge same-domain faces only after all of the legs have been fused.
  return FuseSweep(aStart, aLegs, true);
}

namespace
{
CylinderGeometry ReadCylinder(const TopoDS_Solid& theCylinder)
{
  if (theCylinder.IsNull())
  {
    throw Standard_ConstructionError("Expected a cylinder");
  }

  TopoDS_Face                   aLateralFace;
  NCollection_List<TopoDS_Face> aCaps;
  for (TopExp_Explorer anExplorer(theCylinder, TopAbs_FACE); anExplorer.More(); anExplorer.Next())
  {
    const TopoDS_Face&        aFace = TopoDS::Face(anExplorer.Current());
    const BRepAdaptor_Surface aSurface(aFace);
    if (aSurface.GetType() == GeomAbs_Cylinder && aLateralFace.IsNull())
    {
      aLateralFace = aFace;
    }
    else if (aSurface.GetType() == GeomAbs_Plane)
    {
      aCaps.Append(aFace);
    }
    else
    {
      throw Standard_ConstructionError("Expected one cylindrical face and two planar caps");
    }
  }
  if (aLateralFace.IsNull() || aCaps.Size() != 2)
  {
    throw Standard_ConstructionError("Expected one cylindrical face and two planar caps");
  }

  const BRepAdaptor_Surface aSurface(aLateralFace);
  const gp_Cylinder         aCylinder = aSurface.Cylinder();
  const double              aHeight   = aSurface.LastVParameter() - aSurface.FirstVParameter();
  if (std::abs(aSurface.LastUParameter() - aSurface.FirstUParameter() - 2.0 * std::acos(-1.0))
        > Precision::Angular()
      || !std::isfinite(aHeight) || aHeight <= Precision::Confusion()
      || aCylinder.Radius() <= Precision::Confusion())
  {
    throw Standard_ConstructionError("Expected a full, finite cylinder");
  }

  const gp_Dir anAxis = aCylinder.Axis().Direction();
  for (const TopoDS_Face& aCap : aCaps)
  {
    const gp_Dir aNormal = BRepAdaptor_Surface(aCap).Plane().Axis().Direction();
    if (!aNormal.IsParallel(anAxis, Precision::Angular()))
    {
      throw Standard_ConstructionError("Expected caps perpendicular to the cylinder axis");
    }
  }
  const gp_Pnt aBase = aCylinder.Location().Translated(aSurface.FirstVParameter() * gp_Vec(anAxis));
  const gp_Pnt aFirstCapOrigin = BRepAdaptor_Surface(aCaps.First()).Plane().Location();
  const bool   isFirstCapLower = gp_Vec(aBase, aFirstCapOrigin).Dot(gp_Vec(anAxis)) < 0.5 * aHeight;
  return {aBase,
          anAxis,
          aCylinder.Radius(),
          aHeight,
          isFirstCapLower ? aCaps.First() : aCaps.Last(),
          isFirstCapLower ? aCaps.Last() : aCaps.First()};
}

TopoDS_Solid SweepCylinderTransversely(const TopoDS_Solid& theCylinder,
                                       const gp_Vec&       theTranslation)
{
  const CylinderGeometry aCylinder = ReadCylinder(theCylinder);
  const gp_Dir&          anAxis    = aCylinder.Axis;
  const gp_Vec           anAxisVector(anAxis);
  const double           anAxialMotion = theTranslation.Dot(anAxisVector);
  if (std::abs(anAxialMotion) > Precision::Confusion())
  {
    throw Standard_ConstructionError("Expected translation perpendicular to the cylinder axis");
  }
  const gp_Vec aMotion = theTranslation - anAxialMotion * anAxisVector;
  const double aTravel = aMotion.Magnitude();
  if (aTravel <= Precision::Confusion())
  {
    return theCylinder;
  }

  // The diameter perpendicular to the motion sweeps the stadium's rectangular middle.
  const gp_Pnt& aBase = aCylinder.Base;
  gp_Ax2        aFrame(aBase, anAxis, gp_Dir(aMotion));
  const double  aRadius = aCylinder.Radius;
  aFrame.SetLocation(aBase.Translated(-aRadius * gp_Vec(aFrame.YDirection())));

  NCollection_List<TopoDS_Shape> aTools;
  aTools.Append(BRepPrimAPI_MakeBox(aFrame, aTravel, 2.0 * aRadius, aCylinder.Height).Solid());
  aTools.Append(SolidAtOffset(theCylinder, aMotion));
  return FuseSweep(theCylinder, aTools);
}
} // namespace

TopoDS_Solid BRepAlgoAPI_SolidSweep::Sweep(const TopoDS_Solid&        theSolid,
                                           const CircularTranslation& theMotion)
{
  const gp_Ax1& aPathAxis = theMotion.Axis;
  const double  anAngle   = theMotion.Angle;
  if (!std::isfinite(anAngle) || std::abs(anAngle) > std::acos(-1.0))
  {
    throw Standard_ConstructionError("Expected an arc angle within [-pi, pi]");
  }
  const CylinderGeometry aCylinder      = ReadCylinder(theSolid);
  const gp_Dir&          aPathDirection = aPathAxis.Direction();
  if (!aCylinder.Axis.IsParallel(aPathDirection, Precision::Angular()))
  {
    throw Standard_ConstructionError("Expected parallel path and cylinder axes");
  }

  const gp_Vec aPathDirectionVector(aPathDirection);
  const gp_Vec aFromAxis(aPathAxis.Location(), aCylinder.Base);
  const gp_Vec aRadial     = aFromAxis - aFromAxis.Dot(aPathDirectionVector) * aPathDirectionVector;
  const double aPathRadius = aRadial.Magnitude();
  if (!std::isfinite(aPathRadius))
  {
    throw Standard_ConstructionError("Expected a finite path radius");
  }
  if (std::abs(anAngle) * aPathRadius <= Precision::Confusion())
  {
    return theSolid;
  }

  // Only the part of the radial diameter on the same ray contributes to the
  // sector. The endpoint cylinders cover the footprint outside the arc angles.
  const gp_Vec aHalfWidth = (aCylinder.Radius / aPathRadius) * aRadial;
  const gp_Vec aRise      = aCylinder.Height * gp_Vec(aCylinder.Axis);
  const gp_Pnt anInner    = aCylinder.Base.Translated(
    aPathRadius - aCylinder.Radius <= Precision::Confusion() ? -aRadial : -aHalfWidth);
  const gp_Pnt               anOuter = aCylinder.Base.Translated(aHalfWidth);
  BRepBuilderAPI_MakePolygon aContour;
  aContour.Add(anInner);
  aContour.Add(anOuter);
  aContour.Add(anOuter.Translated(aRise));
  aContour.Add(anInner.Translated(aRise));
  aContour.Close();
  const TopoDS_Face aSection = BRepBuilderAPI_MakeFace(aContour.Wire()).Face();

  gp_Ax1 aRevolutionAxis = aPathAxis;
  if (anAngle < 0.0)
  {
    aRevolutionAxis.Reverse();
  }
  BRepPrimAPI_MakeRevol aSector(aSection, aRevolutionAxis, std::abs(anAngle), true);
  const gp_Pnt          anEnd = aCylinder.Base.Rotated(aPathAxis, anAngle);

  NCollection_List<TopoDS_Shape> aTools;
  aTools.Append(aSector.Shape());
  aTools.Append(SolidAtOffset(theSolid, gp_Vec(aCylinder.Base, anEnd)));
  return FuseSweep(theSolid, aTools);
}

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
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepLib.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <Geom_BSplineSurface.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <NCollection_List.hxx>
#include <Precision.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <Standard_ConstructionError.hxx>
#include <Standard_Failure.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Pln.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

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

namespace
{
struct HelixFrame
{
  gp_Pnt                base;
  gp_Vec                velocity;
  std::array<gp_Vec, 5> normal;
  std::array<gp_Vec, 5> derivative;
};

struct HelixControl
{
  gp_Pnt point;
  gp_Vec derivative;
};

// The two silhouette lines and the two contributing cap semicircles use
// the same Hermite interpolation along the helix, so adjacent patches share
// the same fitted edge. Each semicircle is represented by two rational quarters.
struct HelixSpan
{
  const CylinderGeometry& Cylinder;
  gp_Vec                  Radial;
  gp_Vec                  Tangent;
  gp_Vec                  Advance;
  double                  Begin;
  double                  End;
  bool                    PositiveAxial;
  int                     Intervals = 0;
  double                  Step      = 0.0;

  HelixFrame                       frameAt(double t) const;
  HelixControl                     controlAt(const HelixFrame& frame, int patch, int index) const;
  occ::handle<Geom_BSplineSurface> makePatch(int patch) const;
  gp_Pnt       analyticPoint(const HelixFrame& frame, int patch, double u) const;
  void         addEndpoint(BRepBuilderAPI_Sewing& sewing, const HelixFrame& frame, bool last) const;
  TopoDS_Solid build(double tolerance);
};

HelixFrame HelixSpan::frameAt(double t) const
{
  const double angle        = Begin + t;
  const gp_Vec radial       = std::cos(angle) * Radial + std::sin(angle) * Tangent;
  const gp_Pnt base         = Cylinder.Base.Translated(radial - Radial + angle * Advance);
  const gp_Vec velocity     = -std::sin(angle) * Radial + std::cos(angle) * Tangent + Advance;
  const gp_Vec acceleration = -radial;
  const gp_Vec axis(Cylinder.Axis);
  const gp_Vec transverse = velocity - velocity.Dot(axis) * axis;
  if (transverse.Magnitude() <= Precision::Confusion())
  {
    throw Standard_ConstructionError("Helical sweep has an undefined lateral silhouette");
  }
  const gp_Vec transverseDerivative = acceleration - acceleration.Dot(axis) * axis;
  const gp_Vec u                    = transverse.Normalized();
  const gp_Vec du =
    (transverseDerivative - u.Dot(transverseDerivative) * u) / transverse.Magnitude();
  const gp_Vec w  = axis.Crossed(u);
  const gp_Vec dw = axis.Crossed(du);
  return {base, velocity, {u, w, -u, -w, u}, {du, dw, -du, -dw, du}};
}

HelixControl HelixSpan::controlAt(const HelixFrame& frame, int patch, int index) const
{
  gp_Vec direction;
  gp_Vec derivative;
  double elevation;
  if (patch < 4)
  {
    direction  = frame.normal[patch + (index == 3 ? 1 : 0)];
    derivative = frame.derivative[patch + (index == 3 ? 1 : 0)];
    if (index == 2)
    {
      direction += frame.normal[patch + 1];
      derivative += frame.derivative[patch + 1];
    }
    elevation = ((patch == 0 || patch == 3) == PositiveAxial) ? 0.0 : Cylinder.Height;
  }
  else
  {
    const int seam = patch == 4 ? 1 : 3;
    direction      = frame.normal[seam];
    derivative     = frame.derivative[seam];
    elevation      = index == 1 ? 0.0 : Cylinder.Height;
  }
  return {frame.base.Translated(Cylinder.Radius * direction + elevation * gp_Vec(Cylinder.Axis)),
          frame.velocity + Cylinder.Radius * derivative};
}

occ::handle<Geom_BSplineSurface> HelixSpan::makePatch(int patch) const
{
  const int                  degree = patch < 4 ? 2 : 1;
  NCollection_Array2<gp_Pnt> poles(1, degree + 1, 1, 3 * Intervals + 1);
  NCollection_Array2<double> weights(1, degree + 1, 1, 3 * Intervals + 1);
  NCollection_Array1<double> uk(1, 2), vk(1, Intervals + 1);
  NCollection_Array1<int>    um(1, 2), vm(1, Intervals + 1);
  uk(1) = 0.0;
  uk(2) = 1.0;
  um.Init(degree + 1);
  for (int i = 0; i <= Intervals; ++i)
  {
    vk(i + 1) = i * Step;
    vm(i + 1) = (i == 0 || i == Intervals) ? 4 : 3;
  }
  for (int leg = 0; leg < Intervals; ++leg)
  {
    const HelixFrame start = frameAt(leg * Step);
    const HelixFrame end   = frameAt((leg + 1) * Step);
    for (int i = 1; i <= degree + 1; ++i)
    {
      const HelixControl first  = controlAt(start, patch, i);
      const HelixControl last   = controlAt(end, patch, i);
      const int          column = 3 * leg + 1;
      poles(i, column)          = first.point;
      poles(i, column + 1)      = first.point.Translated((Step / 3.0) * first.derivative);
      poles(i, column + 2)      = last.point.Translated((-Step / 3.0) * last.derivative);
      poles(i, column + 3)      = last.point;
      const double weight       = degree == 2 && i == 2 ? std::sqrt(0.5) : 1.0;
      for (int sample = 0; sample < 4; ++sample)
        weights(i, column + sample) = weight;
    }
  }
  return new Geom_BSplineSurface(poles, weights, uk, vk, um, vm, degree, 3);
}

gp_Pnt HelixSpan::analyticPoint(const HelixFrame& frame, int patch, double u) const
{
  if (patch >= 4)
  {
    const HelixControl first = controlAt(frame, patch, 1);
    return first.point.Translated((u * Cylinder.Height) * gp_Vec(Cylinder.Axis));
  }
  const double b0 = (1.0 - u) * (1.0 - u);
  const double b1 = std::sqrt(2.0) * u * (1.0 - u);
  const double b2 = u * u;
  return gp_Pnt((b0 * controlAt(frame, patch, 1).point.XYZ()
                 + b1 * controlAt(frame, patch, 2).point.XYZ()
                 + b2 * controlAt(frame, patch, 3).point.XYZ())
                / (b0 + b1 + b2));
}

void HelixSpan::addEndpoint(BRepBuilderAPI_Sewing& sewing, const HelixFrame& frame, bool last) const
{
  const double            pi        = std::acos(-1.0);
  const gp_Dir&           direction = Cylinder.Axis;
  const gp_Dir            x(frame.normal[0]);
  const gp_Cylinder       cylinder(gp_Ax3(frame.base, direction, x), Cylinder.Radius);
  BRepBuilderAPI_MakeFace side(cylinder,
                               last ? -pi / 2.0 : pi / 2.0,
                               last ? pi / 2.0 : 3.0 * pi / 2.0,
                               0.0,
                               Cylinder.Height);
  sewing.Add(side.Face());
  const gp_Pnt capCenter =
    frame.base.Translated(last == PositiveAxial ? Cylinder.Height * gp_Vec(direction) : gp_Vec());
  BRepBuilderAPI_MakeEdge rim(gp_Circ(gp_Ax2(capCenter, direction, x), Cylinder.Radius));
  BRepBuilderAPI_MakeWire outline(rim.Edge());
  BRepBuilderAPI_MakeFace cap(outline.Wire());
  sewing.Add(cap.Face());
}

TopoDS_Solid HelixSpan::build(double tolerance)
{
  const double length = End - Begin;
  Intervals           = std::max(2, static_cast<int>(std::ceil(length / 0.08)));
  std::array<occ::handle<Geom_BSplineSurface>, 6> surfaces;
  for (;;)
  {
    Step             = length / Intervals;
    double deviation = 0.0;
    for (int patch = 0; patch < 6; ++patch)
    {
      surfaces[patch] = makePatch(patch);
      for (int interval = 0; interval < Intervals; ++interval)
      {
        for (const double fraction : {0.25, 0.5, 0.75})
        {
          const double     t     = (interval + fraction) * Step;
          const HelixFrame frame = frameAt(t);
          for (int i = 0; i <= 8; ++i)
          {
            const double u = i / 8.0;
            deviation =
              std::max(deviation,
                       analyticPoint(frame, patch, u).Distance(surfaces[patch]->Value(u, t)));
          }
        }
      }
    }
    if (deviation <= tolerance)
    {
      break;
    }
    if (Intervals >= 1024)
    {
      throw Standard_ConstructionError("Helical envelope did not reach the fitting tolerance");
    }
    Intervals *= 2;
  }

  BRepBuilderAPI_Sewing sewing(10.0 * Precision::Confusion());
  for (const auto& surface : surfaces)
  {
    BRepBuilderAPI_MakeFace face(surface, Precision::Confusion());
    if (!face.IsDone())
    {
      throw Standard_ConstructionError("Unable to build a helical envelope face");
    }
    sewing.Add(face.Face());
  }
  addEndpoint(sewing, frameAt(0.0), false);
  addEndpoint(sewing, frameAt(length), true);
  sewing.Perform();
  const TopoDS_Shape& shell = sewing.SewedShape();
  if (shell.IsNull() || shell.ShapeType() != TopAbs_SHELL || sewing.NbFreeEdges() != 0)
  {
    throw Standard_ConstructionError("Helical envelope did not close");
  }
  TopoDS_Solid result = BRepBuilderAPI_MakeSolid(TopoDS::Shell(shell)).Solid();
  if (!BRepLib::OrientClosedSolid(result) || !BRepCheck_Analyzer(result).IsValid())
  {
    throw Standard_ConstructionError("Helical envelope did not produce a valid solid");
  }
  return result;
}
} // namespace

TopoDS_Solid BRepAlgoAPI_SolidSweep::Sweep(const TopoDS_Solid&       theSolid,
                                           const HelicalTranslation& theMotion)
{
  if (!std::isfinite(theMotion.Angle) || !std::isfinite(theMotion.AdvancePerTurn)
      || !std::isfinite(theMotion.SurfaceTolerance) || theMotion.SurfaceTolerance <= 0.0)
  {
    throw Standard_ConstructionError("Expected finite helix parameters and positive fit tolerance");
  }

  const CylinderGeometry cylinder = ReadCylinder(theSolid);
  const double           angle    = std::abs(theMotion.Angle);
  if (angle <= Precision::Angular())
  {
    return theSolid;
  }

  const double pi       = std::acos(-1.0);
  const double fullTurn = 2.0 * pi;
  const gp_Vec pathAxis(theMotion.Axis.Direction());
  const gp_Vec fromAxis(theMotion.Axis.Location(), cylinder.Base);
  const gp_Vec radial  = fromAxis - fromAxis.Dot(pathAxis) * pathAxis;
  const gp_Vec advance = (theMotion.AdvancePerTurn / fullTurn) * pathAxis;
  if (radial.Magnitude() <= Precision::Confusion())
  {
    return Sweep(theSolid, angle * advance);
  }
  const gp_Vec tangent = (theMotion.Angle < 0.0 ? -1.0 : 1.0) * pathAxis.Crossed(radial);
  const gp_Vec axis(cylinder.Axis);

  // The axial velocity is a sinusoid plus a constant. Its sign determines
  // which half of each cap rim contributes to the envelope.
  const double cosine    = tangent.Dot(axis);
  const double sine      = -radial.Dot(axis);
  const double constant  = advance.Dot(axis);
  const double amplitude = std::hypot(cosine, sine);
  if (amplitude <= Precision::Confusion() && std::abs(constant) <= Precision::Confusion())
  {
    throw Standard_ConstructionError("Helical envelope requires axial motion on each span");
  }

  std::vector<double> boundaries{0.0, angle};
  const double        aPathRadius  = radial.Magnitude();
  const double        anExtraAngle = angle - fullTurn;
  // For parallel axes and a partial second turn, put the two temporary
  // section cuts in the repeated azimuth range if the middle cell can still
  // leave enough angular clearance between its endpoint cylinders.
  const bool useOverlapCuts =
    cylinder.Axis.IsParallel(theMotion.Axis.Direction(), Precision::Angular())
    && std::abs(theMotion.AdvancePerTurn) < cylinder.Height && cylinder.Radius < aPathRadius
    && anExtraAngle > 0.0 && anExtraAngle < pi
    && 0.6 * anExtraAngle > 2.0 * std::asin(cylinder.Radius / aPathRadius);
  if (amplitude > std::abs(constant))
  {
    const double phase = std::atan2(sine, cosine);
    const double root  = std::acos(-constant / amplitude);
    for (const double offset : {phase - root, phase + root})
    {
      for (double t = offset + fullTurn * std::ceil(-offset / fullTurn); t < angle; t += fullTurn)
      {
        if (t > Precision::Angular() && angle - t > Precision::Angular())
        {
          boundaries.push_back(t);
        }
      }
    }
  }
  const auto buildSweep = [&](std::vector<double> theBoundaries, bool theSimplify) {
    std::sort(theBoundaries.begin(), theBoundaries.end());
    TopoDS_Solid                   firstSpan;
    NCollection_List<TopoDS_Shape> otherSpans;
    for (std::size_t i = 1; i < theBoundaries.size(); ++i)
    {
      const double begin = theBoundaries[i - 1];
      const double end   = theBoundaries[i];
      if (end - begin <= Precision::Angular())
      {
        continue;
      }
      const double middle   = 0.5 * (begin + end);
      const bool   positive = cosine * std::cos(middle) + sine * std::sin(middle) + constant > 0.0;
      HelixSpan    span{cylinder, radial, tangent, advance, begin, end, positive};
      const TopoDS_Solid sweptSpan = span.build(theMotion.SurfaceTolerance);
      if (firstSpan.IsNull())
      {
        firstSpan = sweptSpan;
      }
      else
      {
        otherSpans.Append(sweptSpan);
      }
    }
    const TopoDS_Solid result = FuseSweep(firstSpan, otherSpans, theSimplify);
    if (result.IsNull() || !BRepCheck_Analyzer(result).IsValid())
    {
      throw Standard_ConstructionError("Helical sweep union is invalid");
    }
    return result;
  };

  // For parallel axes, each horizontal slice sees an arc shorter than half
  // a turn when the advance exceeds twice the height. With positive inner
  // radius, those slice envelopes do not require overlap trimming.
  if (cylinder.Axis.IsParallel(theMotion.Axis.Direction(), Precision::Angular())
      && std::abs(theMotion.AdvancePerTurn) > 2.0 * cylinder.Height
      && aPathRadius - cylinder.Radius > Precision::Confusion())
  {
    try
    {
      return buildSweep({0.0, angle}, false);
    }
    catch (const Standard_Failure&)
    {
      // Retain the existing split construction if the continuous shell fails.
    }
  }

  if (useOverlapCuts)
  {
    TopoDS_Solid aResult;
    try
    {
      aResult = buildSweep({0.0, 0.8 * anExtraAngle, fullTurn + 0.2 * anExtraAngle, angle}, false);
    }
    catch (const Standard_Failure&)
    {
      // Retry with the regular half-turn cuts if the longer middle cell fails.
    }
    if (!aResult.IsNull())
    {
      try
      {
        ShapeUpgrade_UnifySameDomain aCleanup(aResult);
        aCleanup.SetSafeInputMode(true);
        aCleanup.Build();
        const TopoDS_Shape& aCleaned = aCleanup.Shape();
        if (!aCleaned.IsNull() && aCleaned.ShapeType() == TopAbs_SOLID
            && BRepCheck_Analyzer(aCleaned).IsValid())
        {
          return TopoDS::Solid(aCleaned);
        }
      }
      catch (const Standard_Failure&)
      {
        // Preserve the valid Boolean result if cleanup encounters a bad trim.
      }
      return aResult;
    }
  }
  for (double t = pi; t < angle; t += pi)
  {
    boundaries.push_back(t);
  }
  return buildSweep(boundaries, true);
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

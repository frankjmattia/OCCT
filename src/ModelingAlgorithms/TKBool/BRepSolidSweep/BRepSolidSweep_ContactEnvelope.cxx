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

#include <BRepSolidSweep_Sweep.hxx>
#include "BRepSolidSweep_OrbitIntervals.hxx"
#include "BRepSolidSweep_MovingRibbons.hxx"
#include "BRepSolidSweep_RevolvingMotion.hxx"

#include <BOPAlgo_Builder.hxx>
#include <BOPAlgo_Splitter.hxx>
#include <BOPAlgo_MakerVolume.hxx>
#include <BOPTools_AlgoTools3D.hxx>
#include <IntTools_Context.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepLib.hxx>
#include <BRepTools.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_List.hxx>
#include <TopExp.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS_Shell.hxx>
#include <gp_Trsf.hxx>

namespace
{
constexpr double THE_PI = 3.14159265358979323846;

double normalVelocity(const TopoDS_Face&                     theFace,
                      const BRepAdaptor_Surface&             theSurface,
                      const gp_Pnt2d&                        theUV,
                      const BRepSolidSweep::ContactVelocity& theVelocity)
{
  gp_Pnt aPoint;
  gp_Vec aDu, aDv;
  theSurface.D1(theUV.X(), theUV.Y(), aPoint, aDu, aDv);
  gp_Vec aNormal = aDu.Crossed(aDv);
  if (aNormal.SquareMagnitude() <= gp::Resolution())
    return 0.0;
  aNormal.Normalize();
  if (theFace.Orientation() == TopAbs_REVERSED)
    aNormal.Reverse();
  return aNormal.Dot(theVelocity.At(aPoint));
}

struct ContactGenerator
{
  occ::handle<Geom_BoundedCurve> Curve;
  TopoDS_Edge                    SplitEdge;
  TopoDS_Shape                   Source;
};

std::vector<ContactGenerator> collectGenerators(const TopoDS_Solid&                    theSolid,
                                                const BRepSolidSweep::ContactVelocity& theVelocity,
                                                double                                 theTolerance)
{
  std::vector<ContactGenerator> aGenerators;
  for (TopExp_Explorer anExplorer(theSolid, TopAbs_FACE); anExplorer.More(); anExplorer.Next())
  {
    const auto                aFace = TopoDS::Face(anExplorer.Current());
    const BRepAdaptor_Surface aSurface(aFace);
    if (aSurface.GetType() == GeomAbs_Plane)
    {
      for (const auto& anEdge : BRepSolidSweep::PlanarContact(aFace, theVelocity))
      {
        double     aFirst, aLast;
        const auto aCurve = BRep_Tool::Curve(anEdge, aFirst, aLast);
        aGenerators.push_back({new Geom_TrimmedCurve(aCurve, aFirst, aLast), anEdge, aFace});
      }
    }
    else
    {
      for (const auto& aCurve :
           BRepSolidSweep::FitRevolvedGenerators(aFace, theVelocity, theTolerance))
      {
        BRepBuilderAPI_MakeEdge anEdge(aCurve);
        if (!anEdge.IsDone())
          throw Standard_ConstructionError("Could not bound revolved screw contact");
        aGenerators.push_back({aCurve, anEdge.Edge(), aFace});
      }
    }
  }

  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    anAdjacent;
  TopExp::MapShapesAndUniqueAncestors(theSolid, TopAbs_EDGE, TopAbs_FACE, anAdjacent);
  BRepClass3d_SolidClassifier aStationary(theSolid);
  for (int i = 1; i <= anAdjacent.Extent(); ++i)
  {
    const auto& aFaces = anAdjacent.FindFromIndex(i);
    if (aFaces.Extent() != 2)
      continue;
    const auto anEdge = TopoDS::Edge(anAdjacent.FindKey(i));
    const auto aSpans = BRepSolidSweep::SharpContactSpans(anEdge,
                                                          TopoDS::Face(aFaces.First()),
                                                          TopoDS::Face(aFaces.Last()),
                                                          theVelocity);
    for (const auto& aSpan : aSpans)
    {
      const BRepAdaptor_Curve aContact(anEdge);
      const gp_Pnt            aMiddle   = aContact.Value((aSpan.first + aSpan.second) * 0.5);
      bool                    anExposed = true;
      for (const double dt : {-1e-4, 1e-4})
      {
        const gp_Pnt aRelative = aMiddle.Translated(dt * theVelocity.At(aMiddle));
        aStationary.Perform(aRelative, Precision::Confusion());
        anExposed = anExposed && aStationary.State() == TopAbs_OUT;
      }
      if (!anExposed)
        continue;
      double     aFirst, aLast;
      const auto aCurve = BRep_Tool::Curve(anEdge, aFirst, aLast);
      aGenerators.push_back({new Geom_TrimmedCurve(aCurve, aSpan.first, aSpan.second), {}, anEdge});
    }
  }

  return aGenerators;
}

void addRibbon(BRepBuilderAPI_Sewing& theSewing, const occ::handle<Geom_BSplineSurface>& theSurface)
{
  BRepBuilderAPI_MakeFace aFace(theSurface, Precision::Confusion());
  if (!aFace.IsDone() || !BRepCheck_Analyzer(aFace.Face()).IsValid())
    throw Standard_ConstructionError("Could not build a valid contact ribbon");
  theSewing.Add(aFace.Face());
}

void addCaps(BRepBuilderAPI_Sewing&                 aSewing,
             const TopoDS_Solid&                    theSolid,
             const std::vector<ContactGenerator>&   theGenerators,
             const BRepSolidSweep::ContactVelocity& theVelocity,
             double                                 aTolerance,
             const TopLoc_Location*                 theStart,
             const TopLoc_Location*                 theEnd)
{
  BOPAlgo_Splitter aSplitter;
  aSplitter.SetNonDestructive(true);
  aSplitter.SetFuzzyValue(0.3 * aTolerance);
  for (TopExp_Explorer anExplorer(theSolid, TopAbs_FACE); anExplorer.More(); anExplorer.Next())
    aSplitter.AddArgument(anExplorer.Current());
  for (const auto& aGenerator : theGenerators)
    if (!aGenerator.SplitEdge.IsNull())
      aSplitter.AddTool(aGenerator.SplitEdge);
  aSplitter.Perform();
  if (aSplitter.HasErrors())
    throw Standard_ConstructionError("Could not split endpoint faces at screw contact");

  for (TopExp_Explorer anExplorer(aSplitter.Shape(), TopAbs_FACE); anExplorer.More();
       anExplorer.Next())
  {
    const auto                aFace = TopoDS::Face(anExplorer.Current());
    const BRepAdaptor_Surface aSurface(aFace);
    double                    u0, u1, v0, v1;
    BRepTools::UVBounds(aFace, u0, u1, v0, v1);
    bool          anIncoming = false, anOutgoing = false;
    constexpr int aGrid = 11;
    for (int j = 0; j < aGrid * aGrid; ++j)
    {
      const gp_Pnt2d aUV(u0 + (u1 - u0) * (j % aGrid + 0.5) / aGrid,
                         v0 + (v1 - v0) * (j / aGrid + 0.5) / aGrid);
      if (BRepClass_FaceClassifier(aFace, aUV, Precision::PConfusion()).State() != TopAbs_IN)
        continue;
      const double aSpeed = normalVelocity(aFace, aSurface, aUV, theVelocity);
      anIncoming |= aSpeed < -Precision::Confusion();
      anOutgoing |= aSpeed > Precision::Confusion();
    }
    if (anIncoming && anOutgoing)
      throw Standard_ConstructionError("Endpoint fragment crosses a screw contact");
    if (anIncoming && theStart)
      aSewing.Add(aFace.Moved(*theStart));
    else if (anOutgoing && theEnd)
      aSewing.Add(aFace.Moved(*theEnd));
  }
}

TopoDS_Solid closeEnvelope(BRepBuilderAPI_Sewing& aSewing)
{
  aSewing.Perform();
  const auto& aShell = aSewing.SewedShape();
  if (aShell.IsNull() || aShell.ShapeType() != TopAbs_SHELL || aSewing.NbFreeEdges() != 0
      || aSewing.NbMultipleEdges() != 0)
    throw Standard_ConstructionError("Trimmed screw contact envelope did not close");
  auto aSolid = BRepBuilderAPI_MakeSolid(TopoDS::Shell(aShell)).Solid();
  if (!BRepLib::OrientClosedSolid(aSolid) || !BRepCheck_Analyzer(aSolid).IsValid())
    throw Standard_ConstructionError("Trimmed screw contact envelope is not valid");
  return aSolid;
}

template <class PlacementLaw>
TopoDS_Solid trimMovingEnvelope(BRepBuilderAPI_Sewing& theCandidate,
                                const TopoDS_Solid&    theSolid,
                                const PlacementLaw&    thePlacement,
                                double                 theTolerance)
{
  theCandidate.Perform();
  BOPAlgo_Builder aSplit;
  aSplit.SetNonDestructive(true);
  aSplit.SetFuzzyValue(0.3 * theTolerance);
  aSplit.SetUseOBB(true);
  aSplit.SetRunParallel(true);
  for (TopExp_Explorer anIt(theCandidate.SewedShape(), TopAbs_FACE); anIt.More(); anIt.Next())
    aSplit.AddArgument(anIt.Current());
  aSplit.Perform();
  if (aSplit.HasErrors())
    throw Standard_ConstructionError("Could not intersect moving envelope faces");

  Bnd_Box aBounds;
  BRepBndLib::AddOptimal(theSolid, aBounds, false);
  aBounds.Enlarge(theTolerance);
  BRepClass3d_SolidClassifier aStationary(theSolid);
  const auto                  isInsideAt = [&](const gp_Pnt& thePoint, const gp_Trsf& theInverse) {
    const gp_Pnt aLocal = thePoint.Transformed(theInverse);
    if (aBounds.IsOut(aLocal))
      return false;
    aStationary.Perform(aLocal, theTolerance);
    return aStationary.State() == TopAbs_IN;
  };
  constexpr int        aStationCount = 256;
  std::vector<gp_Trsf> anInverse;
  for (int i = 0; i <= aStationCount; ++i)
    anInverse.push_back(thePlacement(double(i) / aStationCount).Inverted());

  IntCurvesFace_ShapeIntersector anIntersector;
  anIntersector.Load(theSolid, Precision::Confusion());
  const auto isBuriedPoint = [&](const gp_Pnt& thePoint) {
    for (const auto& aPose : anInverse)
      if (isInsideAt(thePoint, aPose))
        return true;

    // A short visit can fall between every sampled pose. Intersect a
    // tolerance-checked fit of the inverse trajectory, then classify the
    // actual placement in every resulting interval.
    for (int aCount = 8; aCount <= 128; aCount *= 2)
    {
      const auto aParameters = BRepSolidSweep::RibbonParameters(aCount);
      occ::handle<NCollection_HArray1<gp_Pnt>> aPoints =
        new NCollection_HArray1<gp_Pnt>(1, aCount + 1);
      for (int i = 0; i <= aCount; ++i)
        aPoints->SetValue(i + 1, thePoint.Transformed(anInverse[i * aStationCount / aCount]));
      if (aPoints->First().Distance(aPoints->Last()) <= gp::Resolution()
          && std::all_of(aPoints->begin(), aPoints->end(), [&](const gp_Pnt& theProbe) {
               return aPoints->First().Distance(theProbe) <= gp::Resolution();
             }))
        return false;
      const auto aCurve     = BRepSolidSweep::InterpolateRibbonCurve(aPoints, aParameters);
      bool       isAccurate = true;
      for (int i = 1; i < 4 * aCount && isAccurate; ++i)
      {
        const double t      = double(i) / (4 * aCount);
        const gp_Pnt aProbe = thePoint.Transformed(thePlacement(t).Inverted());
        isAccurate          = aCurve->Value(t).Distance(aProbe) <= 0.15 * theTolerance;
      }
      if (!isAccurate)
        continue;

      const occ::handle<GeomAdaptor_Curve> anAdaptor = new GeomAdaptor_Curve(aCurve);
      anIntersector.Perform(anAdaptor, 0.0, 1.0);
      if (!anIntersector.IsDone())
        throw Standard_ConstructionError("Could not intersect inverse sweep trajectory");
      std::vector<double> aBreaks{0.0, 1.0};
      for (int i = 1; i <= anIntersector.NbPnt(); ++i)
        aBreaks.push_back(std::clamp(anIntersector.WParameter(i), 0.0, 1.0));
      std::sort(aBreaks.begin(), aBreaks.end());
      for (size_t i = 1; i < aBreaks.size(); ++i)
      {
        const double aMiddle = 0.5 * (aBreaks[i - 1] + aBreaks[i]);
        if (isInsideAt(thePoint, thePlacement(aMiddle).Inverted()))
          return true;
      }
      return false;
    }
    throw Standard_ConstructionError("Inverse sweep trajectory exceeds tolerance");
  };

  // Volumetric cells disambiguate nearly coincident trims. A face with
  // material on both sides must disappear even when the border is tangential.
  BOPAlgo_MakerVolume aVolume;
  aVolume.SetIntersect(false);
  aVolume.SetAvoidInternalShapes(true);
  for (TopExp_Explorer anIt(aSplit.Shape(), TopAbs_FACE); anIt.More(); anIt.Next())
    aVolume.AddArgument(anIt.Current());
  aVolume.Perform();
  if (aVolume.HasErrors())
    throw Standard_ConstructionError("Could not form cells from moving envelope");

  const occ::handle<IntTools_Context> aContext = new IntTools_Context();
  std::vector<TopoDS_Solid>           anOccupied;
  for (TopExp_Explorer anIt(aVolume.Shape(), TopAbs_SOLID); anIt.More(); anIt.Next())
  {
    const auto aCell = TopoDS::Solid(anIt.Current());
    if (!BRepCheck_Analyzer(aCell).IsValid())
      throw Standard_ConstructionError("Moving envelope partition has an invalid cell");
    BRepClass3d_SolidClassifier aClassifier(aCell);
    std::optional<gp_Pnt>       aWitness;
    for (TopExp_Explorer aFaces(aCell, TopAbs_FACE); aFaces.More() && !aWitness; aFaces.Next())
    {
      const auto aFace = TopoDS::Face(aFaces.Current());
      gp_Pnt     aPoint;
      gp_Pnt2d   aUV;
      if (BOPTools_AlgoTools3D::PointInFace(aFace, aPoint, aUV, aContext) != 0)
        continue;
      gp_Vec aDu, aDv;
      BRepAdaptor_Surface(aFace).D1(aUV.X(), aUV.Y(), aPoint, aDu, aDv);
      gp_Vec aNormal = aDu.Crossed(aDv);
      if (aNormal.SquareMagnitude() <= gp::Resolution())
        continue;
      aNormal.Normalize();
      if (aFace.Orientation() == TopAbs_REVERSED)
        aNormal.Reverse();
      for (double aFactor : {8.0, 4.0, 2.0})
      {
        const gp_Pnt aProbe = aPoint.Translated((-aFactor * theTolerance) * aNormal);
        aClassifier.Perform(aProbe, theTolerance);
        if (aClassifier.State() == TopAbs_IN)
        {
          aWitness = aProbe;
          break;
        }
      }
    }
    if (!aWitness)
      throw Standard_ConstructionError("Could not classify an interior point of a sweep cell");
    if (isBuriedPoint(*aWitness))
      anOccupied.push_back(aCell);
  }
  if (anOccupied.empty())
    throw Standard_ConstructionError("Moving envelope has no occupied cells");
  if (anOccupied.size() == 1)
    return anOccupied.front();

  // Cells are built from the same split faces. Their shared boundaries
  // separate occupied material and must cancel before closing the result.
  NCollection_IndexedDataMap<TopoDS_Shape, int, TopTools_ShapeMapHasher> aFaces;
  for (const auto& aCell : anOccupied)
    for (TopExp_Explorer anIt(aCell, TopAbs_FACE); anIt.More(); anIt.Next())
    {
      const auto& aFace   = anIt.Current();
      const int   anIndex = aFaces.FindIndex(aFace);
      if (anIndex == 0)
        aFaces.Add(aFace, 1);
      else
        ++aFaces.ChangeFromIndex(anIndex);
    }
  BRepBuilderAPI_Sewing aSewing(theTolerance);
  for (int i = 1; i <= aFaces.Extent(); ++i)
    if (aFaces.FindFromIndex(i) == 1)
      aSewing.Add(aFaces.FindKey(i));
  return closeEnvelope(aSewing);
}
} // namespace

TopoDS_Solid BRepSolidSweep_BuildContactEnvelope(const TopoDS_Solid&                      theSolid,
                                                 const BRepSolidSweep_Sweep::ScrewMotion& theMotion)
{
  if (theSolid.IsNull() || !BRepCheck_Analyzer(theSolid).IsValid())
    throw Standard_ConstructionError("Expected a valid solid for screw sweep");

  const double          aRise      = theMotion.AdvancePerTurn / (2.0 * THE_PI);
  const double          aTolerance = theMotion.SurfaceTolerance;
  const gp_Vec          anOmega(theMotion.Axis.Direction());
  const auto            aVelocity   = BRepSolidSweep::ContactVelocity::Screw(theMotion.Axis, aRise);
  const auto            aGenerators = collectGenerators(theSolid, aVelocity, 0.2 * aTolerance);
  BRepBuilderAPI_Sewing aSewing(aTolerance);
  for (const auto& aGenerator : aGenerators)
    addRibbon(aSewing,
              BRepSolidSweep::FitScrewContactRibbon(aGenerator.Curve,
                                                    theMotion.Axis,
                                                    aRise,
                                                    theMotion.Angle,
                                                    0.8 * aTolerance));
  gp_Trsf anEnd, anAdvance;
  anEnd.SetRotation(theMotion.Axis, theMotion.Angle);
  anAdvance.SetTranslation((aRise * theMotion.Angle) * anOmega);
  anEnd.PreMultiply(anAdvance);
  const TopLoc_Location aStartLocation, anEndLocation(anEnd);
  addCaps(aSewing, theSolid, aGenerators, aVelocity, aTolerance, &aStartLocation, &anEndLocation);
  return closeEnvelope(aSewing);
}

TopoDS_Solid BRepSolidSweep_Sweep::Sweep(const TopoDS_Solid&    theSolid,
                                         const RevolvingMotion& theParameters)
{
  const BRepSolidSweep::RevolvingMotion aLaw(theParameters.Axis,
                                             theParameters.Drift,
                                             theParameters.InitialRate,
                                             theParameters.RateDenominatorSlope);
  return Sweep(theSolid,
               RigidMotion{[&](double theStation) {
                             const auto aPlacement = aLaw.Placement(theStation);
                             const auto aVelocity  = aLaw.Velocity(theStation);
                             return RigidPose{aPlacement,
                                              aVelocity.Angular,
                                              aVelocity.At(gp::Origin().Transformed(aPlacement))};
                           },
                           0.0,
                           theParameters.Length,
                           theParameters.SurfaceTolerance});
}

TopoDS_Solid BRepSolidSweep_Sweep::Sweep(const TopoDS_Solid& theSolid, const RigidMotion& theMotion)
{
  const double aLength    = theMotion.Last - theMotion.First;
  const double aTolerance = theMotion.SurfaceTolerance;
  if (theSolid.IsNull() || !BRepCheck_Analyzer(theSolid).IsValid() || !theMotion.Evaluate
      || !std::isfinite(aLength) || aLength <= 0.0 || !std::isfinite(aTolerance)
      || aTolerance < Precision::Confusion())
    throw Standard_ConstructionError("Expected a valid solid, motion interval and tolerance");

  const auto poseAt = [&](double theTime) {
    const auto aPose = theMotion.Evaluate(theMotion.First + theTime * aLength);
    if (aPose.Placement.IsNegative()
        || std::abs(aPose.Placement.ScaleFactor() - 1.0) > Precision::Angular()
        || !std::isfinite(aPose.Placement.TranslationPart().SquareModulus()))
      throw Standard_ConstructionError("Sweep motion requires a finite rigid placement");
    return aPose;
  };
  const auto referenceVelocity = [](const RigidPose& thePose) {
    const auto   anInverse = thePose.Placement.Inverted();
    const gp_Vec anAngular = thePose.AngularVelocity.Transformed(anInverse);
    const gp_Vec aLinear   = thePose.OriginVelocity.Transformed(anInverse);
    const double aSpeed =
      anAngular.Magnitude() > Precision::Angular() ? anAngular.Magnitude() : aLinear.Magnitude();
    if (!std::isfinite(aSpeed) || aSpeed <= gp::Resolution()
        || !std::isfinite(anAngular.SquareMagnitude()) || !std::isfinite(aLinear.SquareMagnitude()))
      throw Standard_ConstructionError("Sweep motion requires finite nonzero velocity");
    return BRepSolidSweep::ContactVelocity{gp::Origin(), anAngular / aSpeed, aLinear / aSpeed};
  };
  const auto collect = [&](const BRepSolidSweep::ContactVelocity& theVelocity) {
    return collectGenerators(theSolid, theVelocity, 0.15 * aTolerance);
  };
  const auto                       aStartPose       = poseAt(0.0);
  const auto                       anEndPose        = poseAt(1.0);
  const auto                       aStartVelocity   = referenceVelocity(aStartPose);
  const auto                       anEndVelocity    = referenceVelocity(anEndPose);
  const auto                       aStartGenerators = collect(aStartVelocity);
  const auto                       anEndGenerators  = collect(anEndVelocity);
  std::vector<gp_Pnt>              aJunctionPoints;
  std::vector<std::pair<int, int>> aJunctions;
  const auto                       junction = [&](const gp_Pnt& thePoint) {
    for (size_t i = 0; i < aJunctionPoints.size(); ++i)
      if (aJunctionPoints[i].Distance(thePoint) <= 0.25 * aTolerance)
        return static_cast<int>(i);
    aJunctionPoints.push_back(thePoint);
    return static_cast<int>(aJunctionPoints.size() - 1);
  };
  for (const auto& aGenerator : aStartGenerators)
  {
    const int aFirst = junction(aGenerator.Curve->StartPoint());
    const int aLast  = junction(aGenerator.Curve->EndPoint());
    aJunctions.emplace_back(aFirst, aLast);
  }
  std::vector<int> aDegree(aJunctionPoints.size(), 0);
  for (const auto& [aFirst, aLast] : aJunctions)
  {
    ++aDegree[aFirst];
    ++aDegree[aLast];
  }
  for (int aCount : aDegree)
    if (aCount != 2)
      throw Standard_ConstructionError("Moving contact must form closed loops");

  const auto evaluate = [&](double theTime) {
    const auto aPose = theTime == 0.0 ? aStartPose : theTime == 1.0 ? anEndPose : poseAt(theTime);
    const auto aGenerators = theTime == 0.0   ? aStartGenerators
                             : theTime == 1.0 ? anEndGenerators
                                              : collect(referenceVelocity(aPose));
    if (aGenerators.size() != aJunctions.size())
      throw Standard_ConstructionError("Moving contact topology changed");
    std::vector<BRepSolidSweep::MovingGenerator> aMoved;
    for (size_t i = 0; i < aGenerators.size(); ++i)
    {
      if (!aGenerators[i].Source.IsSame(aStartGenerators[i].Source))
        throw Standard_ConstructionError("Moving contact changed its source");
      const auto aCurve =
        occ::down_cast<Geom_BoundedCurve>(aGenerators[i].Curve->Transformed(aPose.Placement));
      aMoved.push_back({aCurve, aJunctions[i].first, aJunctions[i].second});
    }
    return aMoved;
  };

  BRepBuilderAPI_Sewing aSewing(aTolerance);
  for (const auto& aSurface : BRepSolidSweep::FitMovingRibbons(evaluate, aTolerance))
    addRibbon(aSewing, aSurface);
  const TopLoc_Location aStartLocation(aStartPose.Placement), anEndLocation(anEndPose.Placement);
  addCaps(aSewing,
          theSolid,
          aStartGenerators,
          aStartVelocity,
          aTolerance,
          &aStartLocation,
          nullptr);
  addCaps(aSewing, theSolid, anEndGenerators, anEndVelocity, aTolerance, nullptr, &anEndLocation);
  return trimMovingEnvelope(
    aSewing,
    theSolid,
    [&](double theTime) { return poseAt(theTime).Placement; },
    aTolerance);
}

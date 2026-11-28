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

#include "../BRepSolidSweep/BRepSolidSweep_ChordFit.hxx"
#include <BRepSolidSweep_Sweep.hxx>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepLib.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <Standard_Failure.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Torus.hxx>

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace
{
const double kPi        = std::acos(-1.0);
const double kTolerance = 1.0e-4;

const char* surfaceName(int theType)
{
  switch (static_cast<GeomAbs_SurfaceType>(theType))
  {
    case GeomAbs_Cylinder:
      return "cylinder";
    case GeomAbs_Torus:
      return "torus";
    case GeomAbs_Cone:
      return "cone";
    case GeomAbs_Plane:
      return "plane";
    default:
      return "other";
  }
}

TopoDS_Solid makeRoundedCylinder(const gp_Ax2& theFrame,
                                 double        theCrestRadius,
                                 double        theWidth,
                                 double        theCorner)
{
  const double aRadius = theCrestRadius - theCorner;
  const gp_Vec anAxis(theFrame.Direction());
  const auto aLocation = [&](double theZ) { return theFrame.Location().Translated(theZ * anAxis); };
  const auto anAxes    = [&](double theZ) {
    return gp_Ax3(aLocation(theZ), theFrame.Direction(), theFrame.XDirection());
  };

  BRepBuilderAPI_Sewing aSewing(1.0e-6);
  aSewing.Add(BRepBuilderAPI_MakeFace(gp_Torus(anAxes(theCorner), aRadius, theCorner),
                                      0.0,
                                      2.0 * kPi,
                                      -0.5 * kPi,
                                      0.0)
                .Face());
  aSewing.Add(BRepBuilderAPI_MakeFace(gp_Torus(anAxes(theWidth - theCorner), aRadius, theCorner),
                                      0.0,
                                      2.0 * kPi,
                                      0.0,
                                      0.5 * kPi)
                .Face());
  const double      aLength = theWidth - 2.0 * theCorner;
  const gp_Cylinder aCylinder(anAxes(theCorner), theCrestRadius);
  aSewing.Add(BRepBuilderAPI_MakeFace(aCylinder, 0.0, 2.0 * kPi, 0.0, aLength).Face());

  const auto addCap = [&](double theZ, double theRadius) {
    const gp_Ax2            aFrame(aLocation(theZ), theFrame.Direction(), theFrame.XDirection());
    BRepBuilderAPI_MakeEdge anEdge(gp_Circ(aFrame, theRadius));
    BRepBuilderAPI_MakeWire aWire(anEdge.Edge());
    aSewing.Add(BRepBuilderAPI_MakeFace(aWire.Wire()).Face());
  };
  addCap(0.0, aRadius);
  addCap(theWidth, aRadius);

  aSewing.Perform();
  const TopoDS_Shape& aShell = aSewing.SewedShape();
  if (aShell.IsNull() || aShell.ShapeType() != TopAbs_SHELL || aSewing.NbFreeEdges() != 0)
    return {};
  TopoDS_Solid aResult = BRepBuilderAPI_MakeSolid(TopoDS::Shell(aShell)).Solid();
  if (!BRepLib::OrientClosedSolid(aResult) || !BRepCheck_Analyzer(aResult).IsValid())
    return {};
  return aResult;
}

struct ScheduleGuard
{
  BRepSolidSweep::ScrewStationSchedule Previous;

  explicit ScheduleGuard(BRepSolidSweep::ScrewStationSchedule theNext)
      : Previous(BRepSolidSweep::ActiveScrewStationSchedule())
  {
    BRepSolidSweep::SetScrewStationSchedule(theNext);
  }

  ~ScheduleGuard() { BRepSolidSweep::SetScrewStationSchedule(Previous); }
};

void printFits(const char* theLabel, const std::vector<BRepSolidSweep::ScrewPatchFit>& theFits)
{
  std::cout << theLabel << '\n';
  for (const auto& fit : theFits)
  {
    std::cout << "  " << surfaceName(fit.SurfaceType) << " surface=" << fit.SurfaceMeridian
              << " contact=" << fit.ContactMeridian << " proposed=" << fit.ProposedMeridian << "x"
              << fit.ProposedScrew << " accepted=" << fit.AcceptedMeridian << "x"
              << fit.AcceptedScrew << " attempts=" << fit.Attempts
              << " coarse=" << fit.CoarseDeviation << " dense=" << fit.DenseDeviation
              << " uniformStart=" << fit.UniformMeridian << "x" << fit.UniformScrew << '\n';
  }
}

bool checkFormula()
{
  const double aRadius         = 53.5;
  const double aChordTolerance = 0.005;
  const int aCircle = BRepSolidSweep::ChordIntervals(2.0 * kPi * aRadius, aRadius, aChordTolerance);
  const int aLine =
    BRepSolidSweep::ChordIntervals(12.0, std::numeric_limits<double>::infinity(), kTolerance);
  std::cout << "chord formula circle=" << aCircle << " line=" << aLine << '\n';
  return aCircle == 230 && aLine == 1;
}

std::vector<BRepSolidSweep::ScrewPatchFit> sweep(const TopoDS_Solid&                      theBody,
                                                 const BRepSolidSweep_Sweep::ScrewMotion& theMotion,
                                                 BRepSolidSweep::ScrewStationSchedule theSchedule)
{
  ScheduleGuard      aGuard(theSchedule);
  const TopoDS_Solid aSwept = BRepSolidSweep_Sweep::Sweep(theBody, theMotion);
  if (aSwept.IsNull() || !BRepCheck_Analyzer(aSwept).IsValid())
    throw Standard_Failure("Screw sweep did not return a valid solid");
  return BRepSolidSweep::ScrewPatchFits();
}
} // namespace

int RunScrewChordExperiment()
{
  try
  {
    if (!checkFormula())
    {
      std::cerr << "chord formula did not match the circular sagitta count\n";
      return 1;
    }

    const gp_Ax2       aFrame(gp_Pnt(0.0, 1.0, 0.0), gp_Dir(0.0, 0.0, 1.0), gp_Dir(1.0, 0.0, 0.0));
    const TopoDS_Solid aBody = makeRoundedCylinder(aFrame, 0.5, 0.4, 0.05);
    if (aBody.IsNull())
    {
      std::cerr << "rounded cylinder did not sew\n";
      return 1;
    }

    BRepSolidSweep_Sweep::ScrewMotion aMotion;
    aMotion.Axis             = gp_Ax1(gp_Pnt(), gp_Dir(1.0, 0.0, 0.0));
    aMotion.Angle            = kPi / 6.0;
    aMotion.AdvancePerTurn   = 4.0;
    aMotion.SurfaceTolerance = kTolerance;

    const auto aUniform = sweep(aBody, aMotion, BRepSolidSweep::ScrewStationSchedule::Uniform);
    const auto aCurved  = sweep(aBody, aMotion, BRepSolidSweep::ScrewStationSchedule::Curvature);
    printFits("uniform", aUniform);
    printFits("curvature", aCurved);

    int  aUniformCylinder = 0;
    int  aUniformTorus    = 0;
    int  aCurvedCylinder  = 0;
    int  aCurvedTorus     = 0;
    bool aCheckerGap      = false;
    bool aCylinderContact = true;
    bool aTorusCorner     = true;
    bool aDenseHolds      = true;

    for (const auto& fit : aUniform)
    {
      if (fit.SurfaceType == GeomAbs_Cylinder)
        ++aUniformCylinder;
      if (fit.SurfaceType == GeomAbs_Torus)
        ++aUniformTorus;
      if (fit.CoarseDeviation > kTolerance || fit.DenseDeviation > kTolerance)
        aCheckerGap = true;
      if (fit.Attempts != 1 || fit.AcceptedMeridian != fit.UniformMeridian)
        aCheckerGap = true;
    }
    for (const auto& fit : aCurved)
    {
      const bool aIsCylinder = fit.SurfaceType == GeomAbs_Cylinder;
      const bool aIsTorus    = fit.SurfaceType == GeomAbs_Torus;
      if (aIsCylinder)
        ++aCurvedCylinder;
      if (aIsTorus)
        ++aCurvedTorus;
      if (fit.CoarseDeviation > kTolerance || fit.DenseDeviation > kTolerance)
        aDenseHolds = false;
      if (fit.AcceptedMeridian >= fit.UniformMeridian)
        aDenseHolds = false;
      if (aIsCylinder
          && (fit.SurfaceMeridian != 1 || fit.ContactMeridian <= 2
              || fit.AcceptedMeridian > fit.ContactMeridian * 2))
        aCylinderContact = false;
      if (aIsTorus
          && (fit.SurfaceMeridian <= 1 || fit.AcceptedMeridian > fit.SurfaceMeridian * 2
              || fit.ProposedMeridian != fit.SurfaceMeridian))
        aTorusCorner = false;
    }

    const bool aInitializer = aCurvedCylinder > 0 && aCurvedTorus > 0 && aUniformCylinder > 0
                              && aUniformTorus > 0 && !aCheckerGap && aCylinderContact
                              && aTorusCorner && aDenseHolds;
    if (aCheckerGap)
      std::cout << "decision: the checker is the gap\n";
    else
      std::cout << "decision: uniform mid-knot samples already stayed inside the dense probe\n";
    if (aInitializer)
      std::cout << "decision: curvature is a valid initializer. A torus uses its corner "
                   "radius. A straight face takes its station count from the contact curve.\n";
    else
      std::cout << "decision: curvature initializer did not settle"
                << " cylinder=" << aCylinderContact << " torus=" << aTorusCorner
                << " dense=" << aDenseHolds << " bands=" << aCurvedCylinder << "," << aCurvedTorus
                << '\n';
    return aInitializer ? 0 : 1;
  }
  catch (const Standard_Failure& aFailure)
  {
    std::cerr << "experiment failed: " << aFailure.what() << '\n';
    return 1;
  }
}

#ifdef SCREW_CHORD_STANDALONE
int main()
{
  return RunScrewChordExperiment();
}
#else
  #include <gtest/gtest.h>

TEST(BRepSolidSweep_ChordFit, CurvatureGuessSettlesRoundedCylinder)
{
  EXPECT_EQ(RunScrewChordExperiment(), 0);
}
#endif

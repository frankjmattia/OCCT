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

#include <gtest/gtest.h>

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepOffsetAPI_MakeSolidSweep.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRepLib.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Cone.hxx>
#include <gp_Torus.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Pln.hxx>

#include <algorithm>
#include <cmath>

namespace
{
TopoDS_Solid makeRoundedFrustum(const gp_Ax2& theFrame,
                                double        theCrestRadius,
                                double        theWidth,
                                double        theTaper,
                                double        theBackCorner,
                                double        theFrontCorner)
{
  const double aPi          = std::acos(-1.0);
  const double aBackRadius  = theCrestRadius - theBackCorner;
  const double aFrontRadius = aBackRadius + (theBackCorner - theFrontCorner) / std::cos(theTaper)
                              + (theBackCorner - theWidth + theFrontCorner) * std::tan(theTaper);
  const gp_Vec anAxis(theFrame.Direction());
  const auto aLocation = [&](double theZ) { return theFrame.Location().Translated(theZ * anAxis); };
  const auto anAxes    = [&](double theZ) {
    return gp_Ax3(aLocation(theZ), theFrame.Direction(), theFrame.XDirection());
  };

  BRepBuilderAPI_Sewing aSewing(1.0e-6);
  aSewing.Add(BRepBuilderAPI_MakeFace(gp_Torus(anAxes(theBackCorner), aBackRadius, theBackCorner),
                                      0.0,
                                      2.0 * aPi,
                                      -0.5 * aPi,
                                      theTaper)
                .Face());
  aSewing.Add(BRepBuilderAPI_MakeFace(
                gp_Torus(anAxes(theWidth - theFrontCorner), aFrontRadius, theFrontCorner),
                0.0,
                2.0 * aPi,
                theTaper,
                0.5 * aPi)
                .Face());

  const double aBackZ  = theBackCorner * (1.0 + std::sin(theTaper));
  const double aFrontZ = theWidth - theFrontCorner * (1.0 - std::sin(theTaper));
  const double aLength = (aFrontZ - aBackZ) / std::cos(theTaper);
  if (theTaper == 0.0)
  {
    const gp_Cylinder aCylinder(anAxes(aBackZ), theCrestRadius);
    aSewing.Add(BRepBuilderAPI_MakeFace(aCylinder, 0.0, 2.0 * aPi, 0.0, aLength).Face());
  }
  else
  {
    const gp_Cone aCone(anAxes(aBackZ),
                        -theTaper,
                        aBackRadius + theBackCorner * std::cos(theTaper));
    // Two tangent conical faces exercise a different topology from the cylinder.
    aSewing.Add(BRepBuilderAPI_MakeFace(aCone, 0.0, 2.0 * aPi, 0.0, 0.5 * aLength).Face());
    aSewing.Add(BRepBuilderAPI_MakeFace(aCone, 0.0, 2.0 * aPi, 0.5 * aLength, aLength).Face());
  }

  const auto addCap = [&](double theZ, double theRadius) {
    const gp_Ax2            aFrame(aLocation(theZ), theFrame.Direction(), theFrame.XDirection());
    BRepBuilderAPI_MakeEdge anEdge(gp_Circ(aFrame, theRadius));
    BRepBuilderAPI_MakeWire aWire(anEdge.Edge());
    aSewing.Add(BRepBuilderAPI_MakeFace(aWire.Wire()).Face());
  };
  addCap(0.0, aBackRadius);
  addCap(theWidth, aFrontRadius);

  aSewing.Perform();
  const TopoDS_Shape& aShell = aSewing.SewedShape();
  if (aShell.IsNull() || aShell.ShapeType() != TopAbs_SHELL || aSewing.NbFreeEdges() != 0)
    return {};
  TopoDS_Solid aResult = BRepBuilderAPI_MakeSolid(TopoDS::Shell(aShell)).Solid();
  if (!BRepLib::OrientClosedSolid(aResult))
    return {};
  return aResult;
}
} // namespace

TEST(BRepOffsetAPI_MakeSolidSweepTest, AxisymmetricBodiesAlongScrewMotion)
{
  const double aPi   = std::acos(-1.0);
  const double aTilt = aPi / 7.0;
  const gp_Dir aBodyAxis(std::sin(aTilt), 0.0, std::cos(aTilt));
  const gp_Pnt aBase(170.0, 0.0, 0.0);
  const gp_Ax1 anAxis(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0));
  const double aWidth    = 9.6;
  const double anAdvance = 360.0;
  for (const double aTaper : {0.0, aPi / 9.0})
  {
    SCOPED_TRACE(aTaper);
    const TopoDS_Solid aBody =
      makeRoundedFrustum(gp_Ax2(aBase, aBodyAxis), 39.4, aWidth, aTaper, 0.65, 1.1);
    ASSERT_FALSE(aBody.IsNull());

    const TopoDS_Solid aResult = BRepOffsetAPI_MakeSolidSweep::Sweep(
      aBody,
      BRepOffsetAPI_MakeSolidSweep::ScrewMotion{anAxis, 2.5 * aPi, anAdvance, 0.0001});
    ASSERT_FALSE(aResult.IsNull());
    ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

    const double aStation = 2.25 * aPi;
    const gp_Pnt aCenter =
      aBase.Translated(0.5 * aWidth * gp_Vec(aBodyAxis))
        .Rotated(anAxis, aStation)
        .Translated((anAdvance * aStation / (2.0 * aPi)) * gp_Vec(anAxis.Direction()));
    BRepClass3d_SolidClassifier aClassifier(aResult, aCenter, 1.0e-7);
    EXPECT_EQ(aClassifier.State(), TopAbs_IN);
    aClassifier.Perform(gp_Pnt(0.0, 0.0, aCenter.Z()), 1.0e-7);
    EXPECT_EQ(aClassifier.State(), TopAbs_OUT);
  }
}

TEST(BRepOffsetAPI_MakeSolidSweepTest, BoxAlongX)
{
  const TopoDS_Solid aBody = BRepPrimAPI_MakeBox(2.0, 3.0, 4.0).Solid();

  const TopoDS_Shape aResult = BRepOffsetAPI_MakeSolidSweep::Sweep(aBody, gp_Vec(4.0, 0.0, 0.0));

  ASSERT_EQ(aResult.ShapeType(), TopAbs_SOLID);
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(aResult, aProperties);
  EXPECT_NEAR(aProperties.Mass(), 72.0, 1.0e-8);
}

TEST(BRepOffsetAPI_MakeSolidSweepTest, BoxAlongDiagonal)
{
  const TopoDS_Solid aBody = BRepPrimAPI_MakeBox(2.0, 3.0, 4.0).Solid();

  const TopoDS_Shape aResult = BRepOffsetAPI_MakeSolidSweep::Sweep(aBody, gp_Vec(4.0, -5.0, 0.0));

  ASSERT_EQ(aResult.ShapeType(), TopAbs_SOLID);
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(aResult, aProperties);
  EXPECT_NEAR(aProperties.Mass(), 112.0, 1.0e-8);
}

TEST(BRepOffsetAPI_MakeSolidSweepTest, TetrahedronAlongZ)
{
  const gp_Pnt anOrigin(0.0, 0.0, 0.0);
  const gp_Pnt anX(3.0, 0.0, 0.0);
  const gp_Pnt anY(0.0, 4.0, 0.0);
  const gp_Pnt aZ(0.0, 0.0, 2.0);

  BRepBuilderAPI_Sewing aSewing;
  const auto            addTriangle =
    [&aSewing](const gp_Pnt& theFirst, const gp_Pnt& theSecond, const gp_Pnt& theThird) {
      BRepBuilderAPI_MakePolygon aTriangle(theFirst, theSecond, theThird, true);
      aSewing.Add(BRepBuilderAPI_MakeFace(aTriangle.Wire()).Face());
    };
  addTriangle(anOrigin, anY, anX);
  addTriangle(anOrigin, anX, aZ);
  addTriangle(anOrigin, aZ, anY);
  addTriangle(anX, anY, aZ);
  aSewing.Perform();
  const TopoDS_Shape& aShell = aSewing.SewedShape();
  ASSERT_FALSE(aShell.IsNull());
  ASSERT_EQ(aShell.ShapeType(), TopAbs_SHELL);
  TopoDS_Solid aBody = BRepBuilderAPI_MakeSolid(TopoDS::Shell(aShell)).Solid();
  ASSERT_TRUE(BRepLib::OrientClosedSolid(aBody));
  ASSERT_TRUE(BRepCheck_Analyzer(aBody).IsValid());

  const TopoDS_Solid aResult = BRepOffsetAPI_MakeSolidSweep::Sweep(aBody, gp_Vec(0.0, 0.0, 5.0));
  ASSERT_FALSE(aResult.IsNull());
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  // The tetrahedron has volume 4, and its XY projection has area 6.
  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(aResult, aProperties);
  EXPECT_NEAR(aProperties.Mass(), 4.0 + 5.0 * 6.0, 1.0e-8);
}

TEST(BRepOffsetAPI_MakeSolidSweepTest, BoxAlongRetracingPolyline)
{
  const TopoDS_Solid aBody = BRepPrimAPI_MakeBox(2.0, 3.0, 4.0).Solid();

  NCollection_Sequence<gp_Vec> anOffsets;
  anOffsets.Append(gp_Vec(0.0, 0.0, 0.0));
  anOffsets.Append(gp_Vec(4.0, 0.0, 0.0));
  anOffsets.Append(gp_Vec(1.0, 0.0, 0.0));
  anOffsets.Append(gp_Vec(7.0, 0.0, 0.0));
  const TopoDS_Shape aResult = BRepOffsetAPI_MakeSolidSweep::Sweep(aBody, anOffsets);

  ASSERT_EQ(aResult.ShapeType(), TopAbs_SOLID);
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(aResult, aProperties);
  EXPECT_NEAR(aProperties.Mass(), 108.0, 1.0e-8);
}

TEST(BRepOffsetAPI_MakeSolidSweepTest, BoxAlongCrossingPolyline)
{
  const TopoDS_Solid aBody = BRepPrimAPI_MakeBox(2.0, 3.0, 4.0).Solid();

  NCollection_Sequence<gp_Vec> anOffsets;
  anOffsets.Append(gp_Vec(0.0, 0.0, 0.0));
  anOffsets.Append(gp_Vec(8.0, 0.0, 0.0));
  anOffsets.Append(gp_Vec(8.0, 8.0, 0.0));
  anOffsets.Append(gp_Vec(4.0, 8.0, 0.0));
  anOffsets.Append(gp_Vec(4.0, -4.0, 0.0));
  const TopoDS_Shape aResult = BRepOffsetAPI_MakeSolidSweep::Sweep(aBody, anOffsets);

  ASSERT_EQ(aResult.ShapeType(), TopAbs_SOLID);
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  // The four XY footprints have areas 30, 22, 18 and 30. The three
  // corners and the crossing each share 6 square units, with no triple overlap.
  const double aVolume = (30.0 + 22.0 + 18.0 + 30.0 - 4.0 * 6.0) * 4.0;
  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(aResult, aProperties);
  EXPECT_NEAR(aProperties.Mass(), aVolume, 1.0e-8);
}

TEST(BRepOffsetAPI_MakeSolidSweepTest, CylinderAlongTransverseLine)
{
  const double       aRadius = 2.0;
  const double       aHeight = 3.0;
  const double       aTravel = 5.0;
  const TopoDS_Solid aBody   = BRepPrimAPI_MakeCylinder(aRadius, aHeight).Solid();

  const TopoDS_Solid aResult = BRepOffsetAPI_MakeSolidSweep::Sweep(aBody, gp_Vec(aTravel, 0.0, 0.0));

  ASSERT_EQ(aResult.ShapeType(), TopAbs_SOLID);
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  // The footprint is a disk plus the rectangle between the two semicircular ends.
  const double aVolume = aHeight * (std::acos(-1.0) * aRadius * aRadius + 2.0 * aRadius * aTravel);
  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(aResult, aProperties);
  EXPECT_NEAR(aProperties.Mass(), aVolume, 1.0e-8);
}

TEST(BRepOffsetAPI_MakeSolidSweepTest, CylinderAlongAxis)
{
  const double       aRadius = 2.0;
  const double       aHeight = 3.0;
  const double       aTravel = 5.0;
  const TopoDS_Solid aBody   = BRepPrimAPI_MakeCylinder(aRadius, aHeight).Solid();
  const TopoDS_Solid aResult = BRepOffsetAPI_MakeSolidSweep::Sweep(aBody, gp_Vec(0.0, 0.0, aTravel));

  ASSERT_FALSE(aResult.IsNull());
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(aResult, aProperties);
  EXPECT_NEAR(aProperties.Mass(),
              std::acos(-1.0) * aRadius * aRadius * (aHeight + aTravel),
              1.0e-8);
}

TEST(BRepOffsetAPI_MakeSolidSweepTest, CylinderAlongObliqueLine)
{
  const TopoDS_Solid aBody   = BRepPrimAPI_MakeCylinder(2.0, 3.0).Solid();
  const TopoDS_Solid aResult = BRepOffsetAPI_MakeSolidSweep::Sweep(aBody, gp_Vec(5.0, 0.0, 4.0));
  ASSERT_FALSE(aResult.IsNull());
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(aResult, aProperties);
  EXPECT_NEAR(aProperties.Mass(), 28.0 * std::acos(-1.0) + 60.0, 1.0e-8);

  // A horizontal leg would fill this point before the cylinder rises.
  BRepClass3d_SolidClassifier aClassifier(aResult, gp_Pnt(6.0, 0.0, 1.0), 1.0e-7);
  EXPECT_EQ(aClassifier.State(), TopAbs_OUT);
}

TEST(BRepOffsetAPI_MakeSolidSweepTest, CylinderAlongObliqueRetrace)
{
  const TopoDS_Solid           aBody = BRepPrimAPI_MakeCylinder(2.0, 3.0).Solid();
  const gp_Vec                 anEnd(10.0, 0.0, 8.0);
  NCollection_Sequence<gp_Vec> anOffsets;
  anOffsets.Append(gp_Vec(0.0, 0.0, 0.0));
  anOffsets.Append(0.5 * anEnd);
  anOffsets.Append(0.25 * anEnd);
  anOffsets.Append(anEnd);

  const TopoDS_Solid aResult    = BRepOffsetAPI_MakeSolidSweep::Sweep(aBody, anOffsets);
  const TopoDS_Solid aReference = BRepOffsetAPI_MakeSolidSweep::Sweep(aBody, anEnd);
  ASSERT_FALSE(aResult.IsNull());
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(aResult, aProperties);
  EXPECT_NEAR(aProperties.Mass(), 44.0 * std::acos(-1.0) + 120.0, 1.0e-8);

  // Each cut must leave zero material. Equal volume alone cannot detect a shift.
  BRepAlgoAPI_Cut anExtra(aResult, aReference);
  BRepAlgoAPI_Cut aMissing(aReference, aResult);
  ASSERT_FALSE(anExtra.HasErrors());
  ASSERT_FALSE(aMissing.HasErrors());
  EXPECT_FALSE(TopExp_Explorer(anExtra.Shape(), TopAbs_SOLID).More());
  EXPECT_FALSE(TopExp_Explorer(aMissing.Shape(), TopAbs_SOLID).More());
}

TEST(BRepOffsetAPI_MakeSolidSweepTest, CylinderAlongCrossingPolyline)
{
  const TopoDS_Solid           aBody = BRepPrimAPI_MakeCylinder(1.0, 3.0).Solid();
  NCollection_Sequence<gp_Vec> anOffsets;
  anOffsets.Append(gp_Vec(0.0, 0.0, 0.0));
  anOffsets.Append(gp_Vec(8.0, 0.0, 0.0));
  anOffsets.Append(gp_Vec(8.0, 8.0, 0.0));
  anOffsets.Append(gp_Vec(4.0, 8.0, 0.0));
  anOffsets.Append(gp_Vec(4.0, -4.0, 0.0));
  const TopoDS_Solid aResult = BRepOffsetAPI_MakeSolidSweep::Sweep(aBody, anOffsets);
  ASSERT_FALSE(aResult.IsNull());
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  // The crossing is interior, all four legs survive, and the loop stays open inside.
  BRepClass3d_SolidClassifier aClassifier(aResult);
  for (const gp_Pnt& aPoint : {gp_Pnt(4.0, 0.0, 1.5),
                               gp_Pnt(2.0, 0.0, 1.5),
                               gp_Pnt(8.0, 4.0, 1.5),
                               gp_Pnt(6.0, 8.0, 1.5),
                               gp_Pnt(4.0, -3.0, 1.5)})
  {
    aClassifier.Perform(aPoint, 1.0e-7);
    EXPECT_EQ(aClassifier.State(), TopAbs_IN);
  }
  aClassifier.Perform(gp_Pnt(6.0, 4.0, 1.5), 1.0e-7);
  EXPECT_EQ(aClassifier.State(), TopAbs_OUT);

  int          aHorizontalFaces = 0;
  const gp_Dir aVertical(0.0, 0.0, 1.0);
  for (TopExp_Explorer anExplorer(aResult, TopAbs_FACE); anExplorer.More(); anExplorer.Next())
  {
    const BRepAdaptor_Surface aSurface(TopoDS::Face(anExplorer.Current()));
    if (aSurface.GetType() == GeomAbs_Plane
        && aSurface.Plane().Axis().Direction().IsParallel(aVertical, 1.0e-12))
    {
      ++aHorizontalFaces;
    }
  }
  EXPECT_EQ(aHorizontalFaces, 2);
}

TEST(BRepOffsetAPI_MakeSolidSweepTest, TiltedCylinderAlongHelix)
{
  const double       aPi   = std::acos(-1.0);
  const double       aTilt = 25.0 * aPi / 180.0;
  const gp_Dir       aBodyAxis(std::sin(aTilt), 0.0, std::cos(aTilt));
  const gp_Ax1       aPathAxis(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0));
  const TopoDS_Solid aBody =
    BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(8.0, 0.0, 0.0), aBodyAxis), 2.0, 3.0).Solid();

  const TopoDS_Solid aHalfTurn =
    BRepOffsetAPI_MakeSolidSweep::Sweep(aBody,
                                  BRepOffsetAPI_MakeSolidSweep::HelicalTranslation{aPathAxis, aPi, 16.0});
  ASSERT_FALSE(aHalfTurn.IsNull());
  ASSERT_TRUE(BRepCheck_Analyzer(aHalfTurn).IsValid());
  GProp_GProps aHalfTurnProperties;
  BRepGProp::VolumeProperties(aHalfTurn, aHalfTurnProperties);
  EXPECT_NEAR(aHalfTurnProperties.Mass(), 386.521, 0.01); // Standalone prototype baseline.

  const TopoDS_Solid aLongSweep = BRepOffsetAPI_MakeSolidSweep::Sweep(
    aBody,
    BRepOffsetAPI_MakeSolidSweep::HelicalTranslation{aPathAxis, 2.5 * aPi, 16.0});
  ASSERT_FALSE(aLongSweep.IsNull());
  ASSERT_TRUE(BRepCheck_Analyzer(aLongSweep).IsValid());
  GProp_GProps aLongSweepProperties;
  BRepGProp::VolumeProperties(aLongSweep, aLongSweepProperties);
  EXPECT_GT(aLongSweepProperties.Mass(), 2.0 * aHalfTurnProperties.Mass());

  // The center at 405 degrees lies beyond the first turn.
  const double                aStation = 2.25 * aPi;
  const gp_Pnt                aBase(8.0 * std::cos(aStation), 8.0 * std::sin(aStation), 18.0);
  const gp_Pnt                aCenter = aBase.Translated(1.5 * gp_Vec(aBodyAxis));
  BRepClass3d_SolidClassifier aClassifier(aLongSweep, aCenter, 1.0e-7);
  EXPECT_EQ(aClassifier.State(), TopAbs_IN);
}

TEST(BRepOffsetAPI_MakeSolidSweepTest, UprightNonOverlappingHelixUsesOneSpan)
{
  const double       aPi = std::acos(-1.0);
  const gp_Ax1       aPathAxis(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0));
  const TopoDS_Solid aBody =
    BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(8.0, 0.0, 0.0), aPathAxis.Direction()), 2.0, 3.0)
      .Solid();
  const TopoDS_Solid aResult = BRepOffsetAPI_MakeSolidSweep::Sweep(
    aBody,
    BRepOffsetAPI_MakeSolidSweep::HelicalTranslation{aPathAxis, 2.5 * aPi, 16.0});
  ASSERT_FALSE(aResult.IsNull());
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  int aFaceCount = 0;
  for (TopExp_Explorer aFace(aResult, TopAbs_FACE); aFace.More(); aFace.Next())
  {
    ++aFaceCount;
  }
  EXPECT_EQ(aFaceCount, 10); // Six envelope patches and four endpoint faces.

  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(aResult, aProperties, 1.0e-9);
  EXPECT_NEAR(aProperties.Mass(), 332.0 * aPi, 0.002);
}

TEST(BRepOffsetAPI_MakeSolidSweepTest, OverlappingHelixHasNoHalfTurnSection)
{
  const double       aPi = std::acos(-1.0);
  const gp_Ax1       aPathAxis(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0));
  const TopoDS_Solid aBody =
    BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(8.0, 0.0, 0.0), aPathAxis.Direction()), 2.0, 3.0)
      .Solid();
  const TopoDS_Solid aResult = BRepOffsetAPI_MakeSolidSweep::Sweep(
    aBody,
    BRepOffsetAPI_MakeSolidSweep::HelicalTranslation{aPathAxis, 2.5 * aPi, 2.0});
  ASSERT_FALSE(aResult.IsNull());
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  bool hasInnerWall = false;
  bool hasOuterWall = false;
  for (TopExp_Explorer aFace(aResult, TopAbs_FACE); aFace.More(); aFace.Next())
  {
    const BRepAdaptor_Surface aSurface(TopoDS::Face(aFace.Current()));
    if (aSurface.GetType() == GeomAbs_Cylinder)
    {
      const double aRadius = aSurface.Cylinder().Radius();
      hasInnerWall         = hasInnerWall || std::abs(aRadius - 6.0) < 1.0e-7;
      hasOuterWall         = hasOuterWall || std::abs(aRadius - 10.0) < 1.0e-7;
    }
  }
  EXPECT_TRUE(hasInnerWall);
  EXPECT_TRUE(hasOuterWall);

  // At the old artificial half-turn cut, both cap rings were on the left.
  for (TopExp_Explorer anEdge(aResult, TopAbs_EDGE); anEdge.More(); anEdge.Next())
  {
    const BRepAdaptor_Curve aCurve(TopoDS::Edge(anEdge.Current()));
    const double            aFirst  = aCurve.FirstParameter();
    const double            aLength = aCurve.LastParameter() - aFirst;
    for (const double aLevel : {1.0, 4.0})
    {
      bool isSection = true;
      for (const double aFraction : {0.25, 0.5, 0.75})
      {
        const gp_Pnt aPoint = aCurve.Value(aFirst + aFraction * aLength);
        isSection = isSection && aPoint.X() < -5.5 && std::abs(aPoint.Z() - aLevel) < 1.0e-6;
      }
      EXPECT_FALSE(isSection);
    }
  }
}

TEST(BRepOffsetAPI_MakeSolidSweepTest, CylinderAlongQuarterCircle)
{
  const double       aRadius     = 2.0;
  const double       aHeight     = 3.0;
  const double       aPathRadius = 8.0;
  const double       aPi         = std::acos(-1.0);
  const double       anAngle     = aPi / 2.0;
  const gp_Ax2       aFrame(gp_Pnt(aPathRadius, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0));
  const TopoDS_Solid aBody = BRepPrimAPI_MakeCylinder(aFrame, aRadius, aHeight).Solid();
  const gp_Ax1       aPathAxis(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0));
  const TopoDS_Shape aResult =
    BRepOffsetAPI_MakeSolidSweep::Sweep(aBody,
                                  BRepOffsetAPI_MakeSolidSweep::CircularTranslation{aPathAxis, anAngle});

  ASSERT_EQ(aResult.ShapeType(), TopAbs_SOLID);
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  // The footprint is an annular sector and two semicircular end caps.
  const double aVolume =
    aHeight * (2.0 * aRadius * aPathRadius * anAngle + aPi * aRadius * aRadius);
  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(aResult, aProperties);
  EXPECT_NEAR(aProperties.Mass(), aVolume, 1.0e-8);
}

TEST(BRepOffsetAPI_MakeSolidSweepTest, CylinderAlongTightSemicircle)
{
  const double aRadius = 2.0;
  const double aHeight = 3.0;
  const double aPi     = std::acos(-1.0);
  const gp_Ax1 aPathAxis(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0));

  for (const double aPathRadius : {1.0, 2.0})
  {
    SCOPED_TRACE(aPathRadius);
    const gp_Ax2       aFrame(gp_Pnt(aPathRadius, 0.0, 0.0), aPathAxis.Direction());
    const TopoDS_Solid aBody = BRepPrimAPI_MakeCylinder(aFrame, aRadius, aHeight).Solid();
    const TopoDS_Shape aResult =
      BRepOffsetAPI_MakeSolidSweep::Sweep(aBody,
                                    BRepOffsetAPI_MakeSolidSweep::CircularTranslation{aPathAxis, aPi});

    ASSERT_EQ(aResult.ShapeType(), TopAbs_SOLID);
    ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

    // Above the X axis the footprint is a half disk. Below it, it is
    // half the union of the endpoint disks, including their lens overlap.
    const double anOuterRadius = aPathRadius + aRadius;
    const double aDiskOverlap =
      2.0 * aRadius * aRadius * std::acos(aPathRadius / aRadius)
      - 2.0 * aPathRadius * std::sqrt(aRadius * aRadius - aPathRadius * aPathRadius);
    const double anArea =
      0.5 * (aPi * anOuterRadius * anOuterRadius + 2.0 * aPi * aRadius * aRadius - aDiskOverlap);
    GProp_GProps aProperties;
    BRepGProp::VolumeProperties(aResult, aProperties);
    EXPECT_NEAR(aProperties.Mass(), anArea * aHeight, 1.0e-8);
  }
}

namespace
{
void checkRevolvingBox(double theRadialDrift)
{
  const TopoDS_Solid aBox = BRepPrimAPI_MakeBox(gp_Pnt(0.0, 2.0, 0.0), 1.0, 1.0, 1.0).Solid();
  constexpr double aLength = 1.4;
  constexpr double aRate = 0.32;
  constexpr double aSlope = 0.10;
  const BRepOffsetAPI_MakeSolidSweep::RevolvingMotion aMotion{
    gp::OX(), gp_Vec(1.0, theRadialDrift, -0.05), aRate, aSlope, aLength, 1e-5};
  const auto aSweep = BRepOffsetAPI_MakeSolidSweep::Sweep(aBox, aMotion);
  ASSERT_TRUE(BRepCheck_Analyzer(aSweep).IsValid());
  BRepClass3d_SolidClassifier aClassifier(aSweep);
  int anInside = 0, anOutside = 0;
  for (double x : {0.19, 0.6, 1.2, 1.65, 2.2})
    for (double y : {1.7, 2.26, 2.53, 2.89, 3.26, 3.5})
      for (double z : {-0.2, 0.32, 0.82, 1.27, 1.81, 2.1})
      {
        const gp_Pnt aPoint(x, y, z);
        double aDistance = 10.0;
        for (int i = 0; i <= 1024; ++i)
        {
          const double t = aLength * i / 1024.0;
          const double aPhase = (aRate / aSlope) * std::log(1.0 + aSlope * t);
          const double yLocal = y * std::cos(aPhase) + z * std::sin(aPhase) - theRadialDrift * t;
          const double zLocal = z * std::cos(aPhase) - y * std::sin(aPhase) + 0.05 * t;
          const double xLocal = x - t;
          aDistance = std::min(aDistance, std::max({-xLocal, xLocal - 1.0,
            2.0 - yLocal, yLocal - 3.0, -zLocal, zLocal - 1.0}));
        }
        // Discard samples within the maximum change between adjacent poses.
        if (std::abs(aDistance) < 0.003)
          continue;
        aClassifier.Perform(aPoint, 1e-6);
        if (aDistance < 0.0)
        {
          ++anInside;
          EXPECT_EQ(aClassifier.State(), TopAbs_IN) << x << ", " << y << ", " << z;
        }
        else
        {
          ++anOutside;
          EXPECT_EQ(aClassifier.State(), TopAbs_OUT);
        }
      }
  EXPECT_GT(anInside, 20);
  EXPECT_GT(anOutside, 20);
}
} // namespace

TEST(BRepOffsetAPI_MakeSolidSweepTest, RevolvingBoxWithRadialDrift)
{
  checkRevolvingBox(0.50);
}

TEST(BRepOffsetAPI_MakeSolidSweepTest, RevolvingBoxTrimsBuriedContactFold)
{
  checkRevolvingBox(0.18);
}

namespace
{
void checkCurvedRigidSweep(bool theRotates)
{
  const TopoDS_Solid                              aBox = BRepPrimAPI_MakeBox(1.0, 1.0, 1.0).Solid();
  constexpr double                                aLength = 1.4;
  const BRepOffsetAPI_MakeSolidSweep::RigidMotion aMotion{
    [=](double theStation) {
      const double t      = theStation - 2.0;
      const double aPitch = theRotates ? 0.04 * t * t : 0.0;
      const double aYaw   = theRotates ? 0.1 * t : 0.0;
      gp_Trsf      aPitchRotation, aYawRotation, aTranslation;
      aPitchRotation.SetRotation(gp::OY(), aPitch);
      aYawRotation.SetRotation(gp::OZ(), aYaw);
      aTranslation.SetTranslation(
        gp_Vec(1.3 * t - 0.3, 0.45 * t + 0.05 * t * t + 0.16, 0.38 * t + 0.10 * t * t + 0.1));
      const gp_Vec anAngular =
        theRotates ? gp_Vec(-0.08 * t * std::sin(aYaw), 0.08 * t * std::cos(aYaw), 0.1) : gp_Vec();
      return BRepOffsetAPI_MakeSolidSweep::RigidPose{aTranslation * aYawRotation * aPitchRotation,
                                                     anAngular,
                                                     gp_Vec(1.3, 0.45 + 0.1 * t, 0.38 + 0.2 * t)};
    },
    2.0,
    2.0 + aLength,
    1e-5};
  const auto aSweep = BRepOffsetAPI_MakeSolidSweep::Sweep(aBox, aMotion);
  ASSERT_TRUE(BRepCheck_Analyzer(aSweep).IsValid());
  for (TopExp_Explorer anIt(aSweep, TopAbs_SHELL); anIt.More(); anIt.Next())
    EXPECT_TRUE(BRep_Tool::IsClosed(anIt.Current()));

  BRepClass3d_SolidClassifier aClassifier(aSweep);
  int                         anInside = 0, anOutside = 0;
  for (double x : {-0.4, 0.12, 0.51, 0.87, 1.27, 1.76, 2.11, 2.3})
    for (double y : {0.17, 0.41, 0.72, 1.13, 1.56, 1.82})
      for (double z : {0.12, 0.37, 0.81, 1.06, 1.42, 1.84})
      {
        double aDistance = 10.0;
        for (int i = 0; i <= 1024; ++i)
        {
          const double t        = aLength * i / 1024.0;
          const double aYaw     = theRotates ? 0.1 * t : 0.0;
          const double aPitch   = theRotates ? 0.04 * t * t : 0.0;
          const double dx       = x - (1.3 * t - 0.3);
          const double dy       = y - (0.45 * t + 0.05 * t * t + 0.16);
          const double dz       = z - (0.38 * t + 0.10 * t * t + 0.1);
          const double aTurnedX = std::cos(aYaw) * dx + std::sin(aYaw) * dy;
          const double aLocalX  = std::cos(aPitch) * aTurnedX - std::sin(aPitch) * dz;
          const double aLocalY  = std::cos(aYaw) * dy - std::sin(aYaw) * dx;
          const double aLocalZ  = std::sin(aPitch) * aTurnedX + std::cos(aPitch) * dz;
          aDistance             = std::min(
            aDistance,
            std::max({-aLocalX, aLocalX - 1.0, -aLocalY, aLocalY - 1.0, -aLocalZ, aLocalZ - 1.0}));
        }
        // Exclude points within the change possible between two sampled poses.
        if (std::abs(aDistance) < 0.004)
          continue;
        aClassifier.Perform(gp_Pnt(x, y, z), 1e-6);
        anInside += aDistance < 0.0;
        anOutside += aDistance > 0.0;
        EXPECT_EQ(aClassifier.State(), aDistance < 0.0 ? TopAbs_IN : TopAbs_OUT)
          << x << ", " << y << ", " << z;
      }
  EXPECT_GT(anInside, 20);
  EXPECT_GT(anOutside, 20);
}
} // namespace

TEST(BRepOffsetAPI_MakeSolidSweepTest, CurvedTranslationWithFixedOrientation)
{
  checkCurvedRigidSweep(false);
}

TEST(BRepOffsetAPI_MakeSolidSweepTest, CurvedMotionWithChangingRotationAxis)
{
  checkCurvedRigidSweep(true);
}

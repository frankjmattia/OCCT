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
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRepLib.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <gp_Ax2.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Pln.hxx>

#include <cmath>

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

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

#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_SolidSweep.hxx>
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
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <gp_Ax2.hxx>
#include <gp_Pln.hxx>

#include <cmath>

TEST(BRepAlgoAPI_SolidSweepTest, BoxAlongX)
{
  const TopoDS_Solid aBody = BRepPrimAPI_MakeBox(2.0, 3.0, 4.0).Solid();

  const TopoDS_Shape aResult = BRepAlgoAPI_SolidSweep::Sweep(aBody, gp_Vec(4.0, 0.0, 0.0));

  ASSERT_EQ(aResult.ShapeType(), TopAbs_SOLID);
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(aResult, aProperties);
  EXPECT_NEAR(aProperties.Mass(), 72.0, 1.0e-8);
}

TEST(BRepAlgoAPI_SolidSweepTest, BoxAlongDiagonal)
{
  const TopoDS_Solid aBody = BRepPrimAPI_MakeBox(2.0, 3.0, 4.0).Solid();

  const TopoDS_Shape aResult = BRepAlgoAPI_SolidSweep::Sweep(aBody, gp_Vec(4.0, -5.0, 0.0));

  ASSERT_EQ(aResult.ShapeType(), TopAbs_SOLID);
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(aResult, aProperties);
  EXPECT_NEAR(aProperties.Mass(), 112.0, 1.0e-8);
}

TEST(BRepAlgoAPI_SolidSweepTest, TetrahedronAlongZ)
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

  const TopoDS_Solid aResult = BRepAlgoAPI_SolidSweep::Sweep(aBody, gp_Vec(0.0, 0.0, 5.0));
  ASSERT_FALSE(aResult.IsNull());
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  // The tetrahedron has volume 4, and its XY projection has area 6.
  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(aResult, aProperties);
  EXPECT_NEAR(aProperties.Mass(), 4.0 + 5.0 * 6.0, 1.0e-8);
}

TEST(BRepAlgoAPI_SolidSweepTest, BoxAlongRetracingPolyline)
{
  const TopoDS_Solid aBody = BRepPrimAPI_MakeBox(2.0, 3.0, 4.0).Solid();

  NCollection_Sequence<gp_Vec> anOffsets;
  anOffsets.Append(gp_Vec(0.0, 0.0, 0.0));
  anOffsets.Append(gp_Vec(4.0, 0.0, 0.0));
  anOffsets.Append(gp_Vec(1.0, 0.0, 0.0));
  anOffsets.Append(gp_Vec(7.0, 0.0, 0.0));
  const TopoDS_Shape aResult = BRepAlgoAPI_SolidSweep::Sweep(aBody, anOffsets);

  ASSERT_EQ(aResult.ShapeType(), TopAbs_SOLID);
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(aResult, aProperties);
  EXPECT_NEAR(aProperties.Mass(), 108.0, 1.0e-8);
}

TEST(BRepAlgoAPI_SolidSweepTest, BoxAlongCrossingPolyline)
{
  const TopoDS_Solid aBody = BRepPrimAPI_MakeBox(2.0, 3.0, 4.0).Solid();

  NCollection_Sequence<gp_Vec> anOffsets;
  anOffsets.Append(gp_Vec(0.0, 0.0, 0.0));
  anOffsets.Append(gp_Vec(8.0, 0.0, 0.0));
  anOffsets.Append(gp_Vec(8.0, 8.0, 0.0));
  anOffsets.Append(gp_Vec(4.0, 8.0, 0.0));
  anOffsets.Append(gp_Vec(4.0, -4.0, 0.0));
  const TopoDS_Shape aResult = BRepAlgoAPI_SolidSweep::Sweep(aBody, anOffsets);

  ASSERT_EQ(aResult.ShapeType(), TopAbs_SOLID);
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  // The four XY footprints have areas 30, 22, 18 and 30. The three
  // corners and the crossing each share 6 square units, with no triple overlap.
  const double aVolume = (30.0 + 22.0 + 18.0 + 30.0 - 4.0 * 6.0) * 4.0;
  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(aResult, aProperties);
  EXPECT_NEAR(aProperties.Mass(), aVolume, 1.0e-8);
}

TEST(BRepAlgoAPI_SolidSweepTest, CylinderAlongTransverseLine)
{
  const double       aRadius = 2.0;
  const double       aHeight = 3.0;
  const double       aTravel = 5.0;
  const TopoDS_Solid aBody   = BRepPrimAPI_MakeCylinder(aRadius, aHeight).Solid();

  const TopoDS_Solid aResult = BRepAlgoAPI_SolidSweep::Sweep(aBody, gp_Vec(aTravel, 0.0, 0.0));

  ASSERT_EQ(aResult.ShapeType(), TopAbs_SOLID);
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  // The footprint is a disk plus the rectangle between the two semicircular ends.
  const double aVolume = aHeight * (std::acos(-1.0) * aRadius * aRadius + 2.0 * aRadius * aTravel);
  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(aResult, aProperties);
  EXPECT_NEAR(aProperties.Mass(), aVolume, 1.0e-8);
}

TEST(BRepAlgoAPI_SolidSweepTest, CylinderAlongAxis)
{
  const double       aRadius = 2.0;
  const double       aHeight = 3.0;
  const double       aTravel = 5.0;
  const TopoDS_Solid aBody   = BRepPrimAPI_MakeCylinder(aRadius, aHeight).Solid();
  const TopoDS_Solid aResult = BRepAlgoAPI_SolidSweep::Sweep(aBody, gp_Vec(0.0, 0.0, aTravel));

  ASSERT_FALSE(aResult.IsNull());
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(aResult, aProperties);
  EXPECT_NEAR(aProperties.Mass(),
              std::acos(-1.0) * aRadius * aRadius * (aHeight + aTravel),
              1.0e-8);
}

TEST(BRepAlgoAPI_SolidSweepTest, CylinderAlongObliqueLine)
{
  const TopoDS_Solid aBody   = BRepPrimAPI_MakeCylinder(2.0, 3.0).Solid();
  const TopoDS_Solid aResult = BRepAlgoAPI_SolidSweep::Sweep(aBody, gp_Vec(5.0, 0.0, 4.0));
  ASSERT_FALSE(aResult.IsNull());
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(aResult, aProperties);
  EXPECT_NEAR(aProperties.Mass(), 28.0 * std::acos(-1.0) + 60.0, 1.0e-8);

  // A horizontal leg would fill this point before the cylinder rises.
  BRepClass3d_SolidClassifier aClassifier(aResult, gp_Pnt(6.0, 0.0, 1.0), 1.0e-7);
  EXPECT_EQ(aClassifier.State(), TopAbs_OUT);
}

TEST(BRepAlgoAPI_SolidSweepTest, CylinderAlongObliqueRetrace)
{
  const TopoDS_Solid           aBody = BRepPrimAPI_MakeCylinder(2.0, 3.0).Solid();
  const gp_Vec                 anEnd(10.0, 0.0, 8.0);
  NCollection_Sequence<gp_Vec> anOffsets;
  anOffsets.Append(gp_Vec(0.0, 0.0, 0.0));
  anOffsets.Append(0.5 * anEnd);
  anOffsets.Append(0.25 * anEnd);
  anOffsets.Append(anEnd);

  const TopoDS_Solid aResult    = BRepAlgoAPI_SolidSweep::Sweep(aBody, anOffsets);
  const TopoDS_Solid aReference = BRepAlgoAPI_SolidSweep::Sweep(aBody, anEnd);
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

TEST(BRepAlgoAPI_SolidSweepTest, CylinderAlongCrossingPolyline)
{
  const TopoDS_Solid           aBody = BRepPrimAPI_MakeCylinder(1.0, 3.0).Solid();
  NCollection_Sequence<gp_Vec> anOffsets;
  anOffsets.Append(gp_Vec(0.0, 0.0, 0.0));
  anOffsets.Append(gp_Vec(8.0, 0.0, 0.0));
  anOffsets.Append(gp_Vec(8.0, 8.0, 0.0));
  anOffsets.Append(gp_Vec(4.0, 8.0, 0.0));
  anOffsets.Append(gp_Vec(4.0, -4.0, 0.0));
  const TopoDS_Solid aResult = BRepAlgoAPI_SolidSweep::Sweep(aBody, anOffsets);
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

TEST(BRepAlgoAPI_SolidSweepTest, CylinderAlongQuarterCircle)
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
    BRepAlgoAPI_SolidSweep::Sweep(aBody,
                                  BRepAlgoAPI_SolidSweep::CircularTranslation{aPathAxis, anAngle});

  ASSERT_EQ(aResult.ShapeType(), TopAbs_SOLID);
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());

  // The footprint is an annular sector and two semicircular end caps.
  const double aVolume =
    aHeight * (2.0 * aRadius * aPathRadius * anAngle + aPi * aRadius * aRadius);
  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(aResult, aProperties);
  EXPECT_NEAR(aProperties.Mass(), aVolume, 1.0e-8);
}

TEST(BRepAlgoAPI_SolidSweepTest, CylinderAlongTightSemicircle)
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
      BRepAlgoAPI_SolidSweep::Sweep(aBody,
                                    BRepAlgoAPI_SolidSweep::CircularTranslation{aPathAxis, aPi});

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

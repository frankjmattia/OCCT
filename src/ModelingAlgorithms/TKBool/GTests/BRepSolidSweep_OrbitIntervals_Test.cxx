// Created for Open CASCADE Technology.
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include "../BRepSolidSweep/BRepSolidSweep_OrbitIntervals.hxx"
#include "../BRepSolidSweep/BRepSolidSweep_RevolvingMotion.hxx"

#include <BRepAdaptor_Curve.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <gp_Torus.hxx>
#include <gp_Cylinder.hxx>
#include <Geom_Circle.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>
#include <gtest/gtest.h>

TEST(BRepSolidSweep_OrbitIntervalsTest, OverlapRetainsExteriorSources)
{
  const auto aResult = BRepSolidSweep::SweepOrbitIntervals({{5.0, 6.0}, {2.0, 2.5}, {1.0, 4.0}},
                                                          1.0);
  ASSERT_EQ(aResult.size(), 1u);
  EXPECT_DOUBLE_EQ(aResult[0].First, 1.0);
  EXPECT_DOUBLE_EQ(aResult[0].Last, 7.0);
  EXPECT_EQ(aResult[0].FirstSource, 2u);
  EXPECT_EQ(aResult[0].LastSource, 0u);
}

TEST(BRepSolidSweep_OrbitIntervalsTest, ReverseSweepPreservesUncoveredGap)
{
  const auto aResult = BRepSolidSweep::SweepOrbitIntervals({{1.0, 2.0}, {8.0, 9.0}}, -2.0);
  ASSERT_EQ(aResult.size(), 2u);
  EXPECT_DOUBLE_EQ(aResult[0].First, -1.0);
  EXPECT_DOUBLE_EQ(aResult[0].Last, 2.0);
  EXPECT_DOUBLE_EQ(aResult[1].First, 6.0);
  EXPECT_DOUBLE_EQ(aResult[1].Last, 9.0);
}

TEST(BRepSolidSweep_OrbitIntervalsTest, PointContactCanBeExposedOrBuried)
{
  const auto aResult =
    BRepSolidSweep::SweepOrbitIntervals({{1.0, 3.0}, {2.0, 2.0}, {9.0, 9.0}}, 4.0);
  ASSERT_EQ(aResult.size(), 2u);
  EXPECT_EQ(aResult[0].FirstSource, 0u);
  EXPECT_EQ(aResult[0].LastSource, 0u);
  EXPECT_DOUBLE_EQ(aResult[1].First, 9.0);
  EXPECT_DOUBLE_EQ(aResult[1].Last, 13.0);
  EXPECT_EQ(aResult[1].FirstSource, 2u);
  EXPECT_EQ(aResult[1].LastSource, 2u);
}


TEST(BRepSolidSweep_OrbitIntervalsTest, ExactHelixIntersectsBox)
{
  const auto aBox = BRepPrimAPI_MakeBox(gp_Pnt(-2.0, 0.5, -4.0), 4.0, 1.0, 8.0).Solid();
  BRepSolidSweep::OrbitIntersector anIntersector(aBox, 1.0e-7);
  const double aPi = std::acos(-1.0);
  const auto aResult = anIntersector.Occupied(gp_Pnt(1.0, 0.0, 0.0),
                                             gp_Ax1(gp_Pnt(), gp_Dir(0, 0, 1)),
                                             0.25, 0.0, 4.0 * aPi);
  ASSERT_EQ(aResult.size(), 2u);
  EXPECT_NEAR(aResult[0].First, aPi / 6.0, 1.0e-7);
  EXPECT_NEAR(aResult[0].Last, 5.0 * aPi / 6.0, 1.0e-7);
  EXPECT_NEAR(aResult[1].First, 2.0 * aPi + aPi / 6.0, 1.0e-7);
  EXPECT_NEAR(aResult[1].Last, 2.0 * aPi + 5.0 * aPi / 6.0, 1.0e-7);
}

TEST(BRepSolidSweep_OrbitIntervalsTest, NegativeRiseOnAxis)
{
  const auto aBox = BRepPrimAPI_MakeBox(gp_Pnt(-1.0, -1.0, -3.0), 2.0, 2.0, 2.0).Solid();
  BRepSolidSweep::OrbitIntersector anIntersector(aBox, 1.0e-7);
  const auto aResult = anIntersector.Occupied(gp_Pnt(),
                                             gp_Ax1(gp_Pnt(), gp_Dir(0, 0, 1)),
                                             -0.5, 0.0, 8.0);
  ASSERT_EQ(aResult.size(), 1u);
  EXPECT_NEAR(aResult[0].First, 2.0, 1.0e-7);
  EXPECT_NEAR(aResult[0].Last, 6.0, 1.0e-7);
}

TEST(BRepSolidSweep_OrbitIntervalsTest, PlanarContactChordSplitsAtLocalFold)
{
  const gp_Ax1 anAxis(gp_Pnt(), gp_Dir(0, 0, 1));
  // On y=0 the x direction is tangent to the chord. With outward normal +X,
  // the zero-velocity chord instead lies on y=0 within an X-normal plane when
  // X is constant. Tilt the plane so the chord can span X while remaining a
  // contact line. Here y=0 and z=-x describe that contact line.
  const gp_Dir aNormal(1, 0, 1);
  const gp_Pnt aFirst(-1, 0, 1);
  const gp_Pnt aLast(3, 0, -3);
  const auto aKept = BRepSolidSweep::ExposedPlanarChord(aFirst, aLast, aNormal, anAxis);
  ASSERT_TRUE(aKept.has_value());
  EXPECT_DOUBLE_EQ(aKept->First, 0.0);
  EXPECT_NEAR(aKept->Last, 0.25, 1e-15);
  const auto aReversed = BRepSolidSweep::ExposedPlanarChord(aLast, aFirst, aNormal, anAxis);
  ASSERT_TRUE(aReversed.has_value());
  EXPECT_NEAR(aReversed->First, 0.75, 1e-15);
  EXPECT_DOUBLE_EQ(aReversed->Last, 1.0);
  EXPECT_FALSE(BRepSolidSweep::ExposedPlanarChord(gp_Pnt(1, 0, -1), aLast, aNormal, anAxis));
}

TEST(BRepSolidSweep_OrbitIntervalsTest, PlanarContactNeutralChordIsRetained)
{
  const gp_Ax1 anAxis(gp_Pnt(), gp_Dir(0, 0, 1));
  const auto aKept = BRepSolidSweep::ExposedPlanarChord(
    gp_Pnt(-1, 0, 0), gp_Pnt(1, 0, 0), gp_Dir(0, 0, 1), anAxis);
  ASSERT_TRUE(aKept.has_value());
  EXPECT_DOUBLE_EQ(aKept->First, 0.0);
  EXPECT_DOUBLE_EQ(aKept->Last, 1.0);
}

TEST(BRepSolidSweep_OrbitIntervalsTest, CurvedContactAccountsForSurfaceAcceleration)
{
  // Cylinder of radius two at (2,0,0), with angular and axial parameters.
  const gp_Dir aNormal(1, 0, 0);
  const gp_Vec aDu(0, 2, 0), aDv(0, 0, 1), aDuu(-2, 0, 0), aZero;
  const gp_Vec aVelocity(0, -2, -3);
  const auto aClearance = [&](double theAcceleration) {
    return BRepSolidSweep::ContactClearanceAcceleration(
      aNormal, aDu, aDv, aDuu, aZero, aZero, aVelocity, gp_Vec(theAcceleration, 0, 0));
  };
  EXPECT_DOUBLE_EQ(aClearance(0), 2);
  EXPECT_DOUBLE_EQ(aClearance(-2), 0);
  EXPECT_DOUBLE_EQ(aClearance(-3), -1);
}

TEST(BRepSolidSweep_OrbitIntervalsTest, OrbitKeyAndLiftPreserveUnwrappedScrewMotion)
{
  const gp_Ax1 anAxis(gp_Pnt(1.0, 2.0, -3.0), gp_Dir(1.0, 2.0, 4.0));
  const gp_Pnt aPoint(3.0, -2.0, 0.1);
  const double anAngle = 19.25;
  for (const double aRise : {-0.7, 0.2})
  {
    const BRepSolidSweep::ScrewOrbitCoordinates aFrame(anAxis, aRise);
    const gp_Pnt aMoved = aPoint.Rotated(anAxis, anAngle)
      .Translated((anAngle * aRise) * gp_Vec(anAxis.Direction()));
    const auto aKey = aFrame.Key(aPoint);
    EXPECT_LT(aFrame.Key(aMoved).Distance(aKey), 1e-12);
    EXPECT_LT(aFrame.Lift(aKey, aFrame.Axial(aMoved)).Distance(aMoved), 1e-12);
    EXPECT_NEAR(aFrame.Lag(aPoint, aMoved), anAngle, 1e-12);
  }
}


TEST(BRepSolidSweep_OrbitIntervalsTest, AnnularCapContactPreservesBore)
{
  const auto across = BRepSolidSweep::AnnularContactSpans(1.0, 2.0, 0.0);
  ASSERT_EQ(across.size(), 2u);
  EXPECT_EQ(across[0], std::make_pair(-2.0, -1.0));
  EXPECT_EQ(across[1], std::make_pair(1.0, 2.0));

  const auto tangent = BRepSolidSweep::AnnularContactSpans(1.0, 2.0, -1.0);
  ASSERT_EQ(tangent.size(), 1u);
  EXPECT_NEAR(tangent[0].first, -std::sqrt(3.0), 1e-15);
  EXPECT_NEAR(tangent[0].second, std::sqrt(3.0), 1e-15);
  EXPECT_TRUE(BRepSolidSweep::AnnularContactSpans(1.0, 2.0, 2.0).empty());
}

TEST(BRepSolidSweep_OrbitIntervalsTest, OffAxisOrbitIntersectsBothSidesOfAnnulus)
{
  const auto outer = BRepPrimAPI_MakeCylinder(2.0, 4.0).Solid();
  const auto inner = BRepPrimAPI_MakeCylinder(1.0, 4.0).Solid();
  BRepAlgoAPI_Cut cut(outer, inner);
  ASSERT_FALSE(cut.HasErrors());
  TopExp_Explorer solids(cut.Shape(), TopAbs_SOLID);
  ASSERT_TRUE(solids.More());
  BRepSolidSweep::OrbitIntersector intersector(TopoDS::Solid(solids.Current()), 1e-7);
  const double turn = 2.0 * std::acos(-1.0);
  const gp_Ax1 axis(gp_Pnt(1.5, 0.0, 0.0), gp_Dir(0, 0, 1));
  const auto intervals = intersector.Occupied(gp_Pnt(2.5, 0.0, 1.0), axis, 0.15, 0.0, turn);
  ASSERT_EQ(intervals.size(), 2u);
  const double outerCrossing = std::acos(0.25);
  const double innerCrossing = std::acos(-0.75);
  EXPECT_NEAR(intervals[0].First, outerCrossing, 1e-7);
  EXPECT_NEAR(intervals[0].Last, innerCrossing, 1e-7);
  EXPECT_NEAR(intervals[1].First, turn - innerCrossing, 1e-7);
  EXPECT_NEAR(intervals[1].Last, turn - outerCrossing, 1e-7);
}


TEST(BRepSolidSweep_OrbitIntervalsTest, PlanarScrewContactClipsBoreAndRadialSector)
{
  const double pi = std::acos(-1.0);
  for (const double angle : {2.0 * pi, pi})
  {
    const auto outer = BRepPrimAPI_MakeCylinder(2.0, 4.0, angle).Solid();
    const auto inner = BRepPrimAPI_MakeCylinder(1.0, 4.0, angle).Solid();
    BRepAlgoAPI_Cut cut(outer, inner);
    ASSERT_FALSE(cut.HasErrors());

    TopoDS_Face cap;
    for (TopExp_Explorer ex(cut.Shape(), TopAbs_FACE); ex.More(); ex.Next())
    {
      const auto face = TopoDS::Face(ex.Current());
      const BRepAdaptor_Surface surface(face);
      if (surface.GetType() == GeomAbs_Plane
          && surface.Plane().Axis().Direction().IsParallel(gp::DZ(), Precision::Angular())
          && surface.Plane().Location().Z() < 1e-7)
      {
        cap = face;
        break;
      }
    }
    ASSERT_FALSE(cap.IsNull());
    // Normal velocity on the cap is y - 0.5 in this screw motion.
    const gp_Ax1 axis(gp_Pnt(0.0, 0.5, 0.0), gp::DX());
    const auto spans = BRepSolidSweep::PlanarScrewContact(cap, axis, 0.25);
    ASSERT_EQ(spans.size(), 2u);
    for (const auto& edge : spans)
    {
      BRepAdaptor_Curve line(edge);
      for (const double t : {line.FirstParameter(), line.LastParameter()})
      {
        const gp_Pnt p = line.Value(t);
        EXPECT_NEAR(p.Y(), 0.5, 1e-7);
        EXPECT_NEAR(p.Z(), 0.0, 1e-7);
        EXPECT_NEAR(std::abs(p.X()), std::abs(p.X()) < 1.0 ? std::sqrt(0.75) : std::sqrt(3.75), 1e-7);
      }
    }
    // The opposite contact line lies outside the half annulus.
    const auto beyondSector = BRepSolidSweep::PlanarScrewContact(
      cap, gp_Ax1(gp_Pnt(0.0, -0.5, 0.0), gp::DX()), 0.25);
    EXPECT_EQ(beyondSector.size(), angle == pi ? 0u : 2u);
  }
}


TEST(BRepSolidSweep_OrbitIntervalsTest, RevolvedContactRespectsCylinderSector)
{
  const double pi = std::acos(-1.0);
  for (const double angle : {2.0 * pi, pi})
  {
    const auto cylinder = BRepPrimAPI_MakeCylinder(2.0, 4.0, angle).Solid();
    TopoDS_Face lateral;
    for (TopExp_Explorer ex(cylinder, TopAbs_FACE); ex.More(); ex.Next())
    {
      const auto face = TopoDS::Face(ex.Current());
      if (BRepAdaptor_Surface(face).GetType() == GeomAbs_Cylinder)
        lateral = face;
    }
    ASSERT_FALSE(lateral.IsNull());

    const auto roots = BRepSolidSweep::RevolvedScrewContactAngles(lateral, gp::OX(), 0.4, 1.2);
    ASSERT_EQ(roots.size(), angle == pi ? 1u : 2u);
    EXPECT_NEAR(roots[0], std::atan2(0.4, 1.2), 1e-13);
    if (angle > pi)
      EXPECT_NEAR(roots[1], roots[0] + pi, 1e-13);
  }
}

TEST(BRepSolidSweep_OrbitIntervalsTest, RevolvedContactSatisfiesNormalVelocityOnTorus)
{
  const double pi = std::acos(-1.0);
  const gp_Torus torus(gp_Ax3(), 3.0, 0.5);
  const auto face = BRepBuilderAPI_MakeFace(torus, 0.0, 2.0 * pi, -0.5 * pi, 0.5 * pi).Face();
  const gp_Ax1 axis(gp_Pnt(0, 0.5, 0.2), gp::DX());
  constexpr double rise = 0.4;
  const gp_Vec omega(axis.Direction());
  BRepAdaptor_Surface surface(face);
  for (const double fraction : {0.35, 0.5, 0.67})
  {
    const double v = surface.FirstVParameter()
      + fraction * (surface.LastVParameter() - surface.FirstVParameter());
    const auto roots = BRepSolidSweep::RevolvedScrewContactAngles(face, axis, rise, v);
    ASSERT_EQ(roots.size(), 2u);
    for (const double u : roots)
    {
      gp_Pnt point;
      gp_Vec du, dv;
      surface.D1(u, v, point, du, dv);
      const auto velocity = omega.Crossed(gp_Vec(axis.Location(), point)) + rise * omega;
      EXPECT_NEAR(du.Crossed(dv).Normalized().Dot(velocity), 0.0, 1e-12);
    }
  }
}

TEST(BRepSolidSweep_OrbitIntervalsTest, SharpContactChangesAtFaceNormalVelocityRoot)
{
  const auto box = BRepPrimAPI_MakeBox(2.0, 2.0, 2.0).Solid();
  TopoDS_Face right, top;
  for (TopExp_Explorer ex(box, TopAbs_FACE); ex.More(); ex.Next())
  {
    const auto face = TopoDS::Face(ex.Current());
    const auto plane = BRepAdaptor_Surface(face).Plane();
    if (plane.Axis().Direction().IsParallel(gp::DX(), Precision::Angular())
        && std::abs(plane.Location().X() - 2.0) < 1e-7)
      right = face;
    if (plane.Axis().Direction().IsParallel(gp::DZ(), Precision::Angular())
        && std::abs(plane.Location().Z() - 2.0) < 1e-7)
      top = face;
  }
  ASSERT_FALSE(right.IsNull());
  ASSERT_FALSE(top.IsNull());
  TopoDS_Edge shared;
  for (TopExp_Explorer a(right, TopAbs_EDGE); a.More(); a.Next())
    for (TopExp_Explorer b(top, TopAbs_EDGE); b.More(); b.Next())
      if (a.Current().IsSame(b.Current()))
        shared = TopoDS::Edge(a.Current());
  ASSERT_FALSE(shared.IsNull());

  // On the right face n.v = rise, while on the top face n.v = y - 1.
  const gp_Ax1 axis(gp_Pnt(0.0, 1.0, 0.0), gp::DX());
  const auto spans = BRepSolidSweep::SharpScrewContactSpans(shared, right, top, axis, 0.25);
  ASSERT_EQ(spans.size(), 1u);
  const BRepAdaptor_Curve section(shared);
  const double a = section.Value(spans[0].first).Y();
  const double b = section.Value(spans[0].second).Y();
  EXPECT_NEAR(std::min(a, b), 0.0, 1e-7);
  EXPECT_NEAR(std::max(a, b), 1.0, 1e-7);
}


TEST(BRepSolidSweep_OrbitIntervalsTest, ScrewRibbonRetainsRationalGeneratorAndMotionTolerance)
{
  const occ::handle<Geom_Circle> circle = new Geom_Circle(gp_Ax2(), 2.0);
  const occ::handle<Geom_TrimmedCurve> generator = new Geom_TrimmedCurve(circle, 0.1, 2.2);
  const gp_Ax1 axis(gp_Pnt(1.1, -0.3, 0.7), gp_Dir(1.0, 2.0, 4.0));
  const gp_Vec direction(axis.Direction());
  constexpr double rise = -0.4, angle = 5.8, tolerance = 1e-6;
  const auto ribbon = BRepSolidSweep::FitScrewContactRibbon(generator, axis, rise, angle, tolerance);
  double u0, u1, v0, v1;
  ribbon->Bounds(u0, u1, v0, v1);
  EXPECT_DOUBLE_EQ(v0, 0.0);
  EXPECT_DOUBLE_EQ(v1, angle);
  ASSERT_TRUE(BRepCheck_Analyzer(BRepBuilderAPI_MakeFace(ribbon, 1e-7).Face()).IsValid());
  for (int i = 0; i <= 14; ++i)
  {
    const double u = u0 + (u1 - u0) * i / 14.0;
    const auto start = ribbon->Value(u, 0.0);
    EXPECT_NEAR(start.XYZ().Modulus(), 2.0, 1e-13);
    for (int j = 0; j <= 201; ++j)
    {
      const double t = angle * j / 201.0;
      const auto expected = start.Rotated(axis, t).Translated((rise * t) * direction);
      EXPECT_LE(ribbon->Value(u, t).Distance(expected), tolerance);
    }
  }
}


TEST(BRepSolidSweep_OrbitIntervalsTest, RevolvedGeneratorStopsAtSectorBoundary)
{
  const gp_Cylinder cylinder(gp_Ax3(), 2.0);
  const auto face = BRepBuilderAPI_MakeFace(cylinder, 0.0, 0.35, 0.0, 4.0).Face();
  const auto curves = BRepSolidSweep::FitRevolvedScrewGenerators(face, gp::OX(), 0.4, 1e-6);
  ASSERT_EQ(curves.size(), 1u);
  const double firstV = 0.4 / std::tan(0.35);
  EXPECT_NEAR(curves[0]->StartPoint().Z(), firstV, 1e-7);
  EXPECT_NEAR(curves[0]->EndPoint().Z(), 4.0, 1e-7);
  for (int i = 0; i <= 76; ++i)
  {
    const auto p = curves[0]->Value(i / 76.0);
    EXPECT_NEAR(std::hypot(p.X(), p.Y()), 2.0, 1e-6);
    EXPECT_NEAR(p.X() * 0.4 - p.Y() * p.Z(), 0.0, 4e-6);
  }
}

TEST(BRepSolidSweep_OrbitIntervalsTest, RigidPlacementScrewReproducesIntermediatePoses)
{
  const gp_Ax1 axis(gp_Pnt(2.0, -4.0, 1.0), gp_Dir(0.4, 0.6, -1.2));
  constexpr double angle = 2.83, advance = -0.5;
  gp_Trsf initial;
  initial.SetTranslation(gp_Vec(1.3, 2.7, -9.0));
  const auto at = [&](double fraction) {
    gp_Trsf rotation, translation;
    rotation.SetRotation(axis, angle * fraction);
    translation.SetTranslation((advance * fraction) * gp_Vec(axis.Direction()));
    return translation * rotation;
  };
  const auto span = BRepSolidSweep::ScrewBetweenPlacements(initial, at(1.0) * initial);
  EXPECT_NEAR(span.Angle, angle, 1e-13);
  EXPECT_NEAR(span.Advance, advance, 1e-13);
  for (const double fraction : {0.0, 0.1, 0.5, 1.0})
    for (const auto& p : {gp_Pnt(), gp_Pnt(1, 2, 3), gp_Pnt(-4, -2, 1)})
      EXPECT_LT(p.Transformed(span.Transformation(fraction)).Distance(p.Transformed(at(fraction))), 1e-12);
}

TEST(BRepSolidSweep_OrbitIntervalsTest, RigidPlacementScrewRetainsPureTranslation)
{
  gp_Trsf start, end;
  start.SetTranslation(gp_Vec(-2, 3, 7));
  end.SetTranslation(gp_Vec(1, -2, 5));
  const auto span = BRepSolidSweep::ScrewBetweenPlacements(start, end);
  EXPECT_DOUBLE_EQ(span.Angle, 0.0);
  EXPECT_LT(gp_Pnt().Transformed(span.Transformation(0.5)).Distance(gp_Pnt(1.5, -2.5, -1.0)), 1e-13);
}


TEST(BRepSolidSweep_OrbitIntervalsTest, RevolvingDriftVelocityMatchesPoseDerivative)
{
  const gp_Ax1 axis(gp_Pnt(1.1, -1.7, 0.3), gp_Dir(1.0, 0.5, -3.0));
  const gp_Pnt p(-2.0, -0.4, 0.7);
  for (const double slope : {-0.3, 0.0, 0.2})
  {
    const BRepSolidSweep::RevolvingMotion motion(axis, gp_Vec(1.0, 0.1, -0.05), 2.31, slope);
    for (const double station : {0.0, 0.4, 1.2})
    {
      constexpr double delta = 1.0e-5;
      const auto before = p.Transformed(motion.Placement(station - delta));
      const auto after = p.Transformed(motion.Placement(station + delta));
      const auto current = p.Transformed(motion.Placement(station));
      const gp_Vec numerical = gp_Vec(before, after) / (2.0 * delta);
      EXPECT_LT((numerical - motion.Velocity(station).At(current)).Magnitude(), 6e-9);
      EXPECT_NEAR(motion.AngularRate(station), 2.31 / (1.0 + slope * station), 1e-14);
    }
  }
}

TEST(BRepSolidSweep_OrbitIntervalsTest, RevolvingDriftWithAxialTranslationReducesToScrew)
{
  const gp_Pnt initial(0.0, 2.0, -1.0);
  const BRepSolidSweep::RevolvingMotion aMotion(gp::OX(), gp_Vec(1, 0, 0), -2.0, 0.0);
  const auto placed = initial.Transformed(aMotion.Placement(0.75));
  const auto expected = initial.Rotated(gp::OX(), -1.5).Translated(gp_Vec(0.75, 0, 0));
  EXPECT_LT(placed.Distance(expected), 1e-13);
}

TEST(BRepSolidSweep_OrbitIntervalsTest, RevolvingDriftInstantaneousScrewPullsBackVelocity)
{
  const gp_Ax1 axis(gp_Pnt(1, 2, 3), gp_Dir(0.5, -1.0, 0.0));
  const gp_Pnt reference(2.0, 1.3, -0.7);
  for (const double rate : {-2.2, 1.7})
  {
    const BRepSolidSweep::RevolvingMotion motion(axis, gp_Vec(1, 0.23, 0.1), rate, 0.1);
    for (const double station : {0.0, 0.5, 1.1})
    {
      const auto screw = motion.ScrewAt(station);
      ASSERT_TRUE(screw.has_value());
      const gp_Vec direction(screw->Axis.Direction());
      const gp_Vec speed = std::abs(motion.AngularRate(station)) *
          (direction.Crossed(gp_Vec(screw->Axis.Location(), reference)) + screw->Rise * direction);
      EXPECT_LT((speed - motion.ReferenceVelocity(station).At(reference)).Magnitude(), 1e-13);
      const auto placement = motion.Placement(station);
      EXPECT_LT((speed.Transformed(placement) -
          motion.Velocity(station).At(reference.Transformed(placement))).Magnitude(), 1e-13);
    }
  }
}

TEST(BRepSolidSweep_OrbitIntervalsTest, TranslationFindsNarrowTorusContactAtTrim)
{
  constexpr double  aMajor = 2.5, aMinor = 0.015;
  const double      aV0 = 1.5 * std::acos(-1.0);
  const TopoDS_Face aFace =
    BRepBuilderAPI_MakeFace(gp_Torus(gp_Ax3(), aMajor, aMinor), -0.18, 0.18, aV0, aV0 + 1.0).Face();
  const BRepSolidSweep::ContactVelocity aVelocity{gp::Origin(), gp_Vec(), gp_Vec(1.0, 0.0, 0.017)};
  const auto aCurves = BRepSolidSweep::FitRevolvedGenerators(aFace, aVelocity, 1e-6);
  ASSERT_EQ(aCurves.size(), 2u);
  EXPECT_LE(aCurves.front()->StartPoint().Distance(aCurves.back()->StartPoint()), 1e-6);
  EXPECT_GT(aCurves.front()->EndPoint().Distance(aCurves.back()->EndPoint()), 0.8);
  for (const auto& aCurve : aCurves)
    for (int i = 0; i <= 32; ++i)
    {
      const auto   aPoint = aCurve->Value(i / 32.0);
      const double u      = std::atan2(aPoint.Y(), aPoint.X());
      const double v      = std::atan2(aPoint.Z(), std::hypot(aPoint.X(), aPoint.Y()) - aMajor);
      EXPECT_NEAR(std::cos(u) * std::cos(v) + 0.017 * std::sin(v), 0.0, 1e-4);
    }
}

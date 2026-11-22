// Created for Open CASCADE Technology.
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#ifndef _BRepSolidSweep_OrbitIntervals_HeaderFile
#define _BRepSolidSweep_OrbitIntervals_HeaderFile

#include <Adaptor3d_CurveOnSurface.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRep_Tool.hxx>
#include <Geom2d_Curve.hxx>
#include <BRepAlgoAPI_Section.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <gp_Pln.hxx>
#include <gp_Trsf.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <Geom2d_Line.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Line.hxx>
#include <Geom_BSplineCurve.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <NCollection_HArray1.hxx>
#include <array>
#include <Geom_BSplineSurface.hxx>
#include <GeomConvert.hxx>
#include <NCollection_Array2.hxx>
#include <IntCurvesFace_ShapeIntersector.hxx>
#include <Precision.hxx>
#include <Standard_ConstructionError.hxx>
#include <TopoDS_Solid.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax3.hxx>
#include <cmath>
#include <algorithm>
#include <cstddef>
#include <vector>
#include <optional>
#include <utility>

namespace BRepSolidSweep
{
//! Principal rigid screw joining two placements. Angles are in [0, pi].
//! For a zero-angle span, only Translation is defined geometrically.
struct ScrewPlacementSpan
{
  gp_Ax1 Axis;
  double Angle;
  double Advance;
  gp_Vec Translation;

  gp_Trsf Transformation(double theFraction) const
  {
    gp_Trsf aResult;
    if (Angle == 0.0)
    {
      aResult.SetTranslation(theFraction * Translation);
      return aResult;
    }

    aResult.SetRotation(Axis, theFraction * Angle);
    gp_Trsf anAdvance;
    anAdvance.SetTranslation((theFraction * Advance) * gp_Vec(Axis.Direction()));
    aResult.PreMultiply(anAdvance);
    return aResult;
  }
};

inline ScrewPlacementSpan ScrewBetweenPlacements(const gp_Trsf& theStart, const gp_Trsf& theEnd)
{
  if (std::abs(theStart.ScaleFactor() - 1.0) > Precision::Confusion()
      || std::abs(theEnd.ScaleFactor() - 1.0) > Precision::Confusion())
    throw Standard_ConstructionError("Screw placements must be rigid and orientation-preserving");

  const gp_Trsf aRelative = theEnd * theStart.Inverted();
  const gp_Vec  aTranslation(aRelative.TranslationPart());
  gp_XYZ        anAxis;
  double        anAngle;
  if (!aRelative.GetRotation(anAxis, anAngle) || std::abs(anAngle) <= Precision::Angular())
    return {gp_Ax1(), 0.0, 0.0, aTranslation};

  if (anAngle < 0.0)
  {
    anAxis.Reverse();
    anAngle = -anAngle;
  }
  const double aPi = std::acos(-1.0);
  if (anAngle > aPi)
  {
    anAxis.Reverse();
    anAngle = 2.0 * aPi - anAngle;
  }

  const gp_Dir aDirection(anAxis);
  const gp_Vec aUnit(aDirection);
  const double anAdvance   = aTranslation.Dot(aUnit);
  const gp_Vec aTransverse = aTranslation - anAdvance * aUnit;
  const gp_Vec aCenter = 0.5 * (aTransverse + aUnit.Crossed(aTransverse) / std::tan(0.5 * anAngle));
  return {gp_Ax1(gp_Pnt(aCenter.XYZ()), aDirection), anAngle, anAdvance, aTranslation};
}

//! Occupied angular interval of the stationary solid on an unwrapped screw orbit.
//! A point contact is represented by equal First and Last values.
struct OrbitInterval
{
  double First;
  double Last;
};

//! Second derivative of signed distance along the inverse trajectory at a
//! regular tangential contact. The surface derivatives and motion derivatives
//! must be expressed in the same frame.
inline double ContactClearanceAcceleration(const gp_Dir& theOutward,
                                           const gp_Vec& theDu,
                                           const gp_Vec& theDv,
                                           const gp_Vec& theDuu,
                                           const gp_Vec& theDuv,
                                           const gp_Vec& theDvv,
                                           const gp_Vec& theVelocity,
                                           const gp_Vec& theAcceleration)
{
  const gp_Vec aCross             = theDu.Crossed(theDv);
  const double aMetricDeterminant = aCross.SquareMagnitude();
  if (!(aMetricDeterminant > gp::Resolution()))
    throw Standard_ConstructionError("Singular surface at screw contact");

  // Recover the parameter velocity of the tangent trajectory using the
  // reciprocal basis. Subtract the surface's normal acceleration.
  const double aU                   = theVelocity.Crossed(theDv).Dot(aCross) / aMetricDeterminant;
  const double aV                   = theDu.Crossed(theVelocity).Dot(aCross) / aMetricDeterminant;
  const gp_Vec aSurfaceAcceleration = aU * aU * theDuu + 2.0 * aU * aV * theDuv + aV * aV * theDvv;
  return gp_Vec(theOutward).Dot(theAcceleration - aSurfaceAcceleration);
}

//! Material spans of a line on an annular face. Coordinates are measured
//! along the line from its closest point to the common center. A tangent
//! to the bore leaves one connected span; a tangent to the outer edge
//! has no positive-length span to sweep.
inline std::vector<std::pair<double, double>> AnnularContactSpans(double theInnerRadius,
                                                                  double theOuterRadius,
                                                                  double theOffset)
{
  if (!std::isfinite(theInnerRadius) || !std::isfinite(theOuterRadius) || !std::isfinite(theOffset)
      || theInnerRadius < 0.0 || theOuterRadius <= theInnerRadius)
    throw Standard_ConstructionError("Invalid annular contact geometry");

  const double aDistance = std::abs(theOffset);
  if (aDistance >= theOuterRadius)
    return {};

  const auto aHalfChord = [&](double theRadius) {
    return std::sqrt((theRadius - aDistance) * (theRadius + aDistance));
  };
  const double anOuter = aHalfChord(theOuterRadius);
  if (aDistance >= theInnerRadius)
    return {{-anOuter, anOuter}};

  const double anInner = aHalfChord(theInnerRadius);
  return {{-anOuter, -anInner}, {anInner, anOuter}};
}

//! Affine velocity in the initial solid's reference frame.
struct ContactVelocity
{
  gp_Pnt Origin;
  gp_Vec Angular;
  gp_Vec Linear;

  gp_Vec At(const gp_Pnt& thePoint) const
  {
    return Angular.Crossed(gp_Vec(Origin, thePoint)) + Linear;
  }

  static ContactVelocity Screw(const gp_Ax1& theAxis, double theRise)
  {
    const gp_Vec anAngular(theAxis.Direction());
    return {theAxis.Location(), anAngular, theRise * anAngular};
  }
};

//! Portions of the zero-normal-velocity line that lie on a trimmed planar face.
//! An identically tangent plane contributes no isolated contact line.
inline std::vector<TopoDS_Edge> PlanarContact(const TopoDS_Face&     theFace,
                                              const ContactVelocity& theVelocity)
{
  const BRepAdaptor_Surface aSurface(theFace);
  if (aSurface.GetType() != GeomAbs_Plane)
    throw Standard_ConstructionError("Expected a planar screw contact face");

  const gp_Pln aPlane = aSurface.Plane();
  const gp_Vec aNormal(aPlane.Axis().Direction());
  const gp_Vec aGradient = aNormal.Crossed(theVelocity.Angular);
  if (aGradient.SquareMagnitude() <= (Precision::Angular() * Precision::Angular()))
    return {};

  const gp_Vec aVelocity     = theVelocity.At(aPlane.Location());
  const gp_Pnt aContactPoint = aPlane.Location().Translated(
    (-aNormal.Dot(aVelocity) / aGradient.SquareMagnitude()) * aGradient);
  BRepAlgoAPI_Section aSection(theFace, gp_Pln(aContactPoint, gp_Dir(aGradient)), false);
  aSection.Build();
  if (!aSection.IsDone() || aSection.HasErrors())
    throw Standard_ConstructionError("Could not intersect planar screw contact with face trims");

  std::vector<TopoDS_Edge> anEdges;
  for (TopExp_Explorer anExplorer(aSection.Shape(), TopAbs_EDGE); anExplorer.More();
       anExplorer.Next())
    anEdges.push_back(TopoDS::Edge(anExplorer.Current()));
  return anEdges;
}

inline std::vector<TopoDS_Edge> PlanarScrewContact(const TopoDS_Face& theFace,
                                                   const gp_Ax1&      theAxis,
                                                   double             theRise)
{
  return PlanarContact(theFace, ContactVelocity::Screw(theAxis, theRise));
}

inline double RevolvedNormalSpeed(const BRepAdaptor_Surface& theSurface,
                                  const ContactVelocity&     theVelocity,
                                  double                     theU,
                                  double                     theV)
{
  gp_Pnt aPoint;
  gp_Vec aDu, aDv;
  theSurface.D1(theU, theV, aPoint, aDu, aDv);
  const gp_Vec aNormal = aDu.Crossed(aDv);
  if (aNormal.SquareMagnitude() <= gp::Resolution())
    throw Standard_ConstructionError("Singular revolved contact station");
  return aNormal.Normalized().Dot(theVelocity.At(aPoint));
}

//! Isolated roots of the normal velocity on a regular cylinder, cone or
//! torus at a given meridian parameter. Roots are expressed in the face's
//! UV chart, including its angular phase and trimming wires.
inline std::array<std::optional<double>, 2> RevolvedContactRoots(const TopoDS_Face&     theFace,
                                                                 const ContactVelocity& theVelocity,
                                                                 double                 theV)
{
  const BRepAdaptor_Surface aSurface(theFace);
  const auto                aType = aSurface.GetType();
  if (aType != GeomAbs_Cylinder && aType != GeomAbs_Cone && aType != GeomAbs_Torus)
    throw Standard_ConstructionError("Expected a cylinder, cone or torus at screw contact");

  const auto   aFirst       = aSurface.FirstUParameter();
  const auto   aLast        = aSurface.LastUParameter();
  const double aPi          = std::acos(-1.0);
  const double aPeriod      = 2.0 * aPi;
  const auto   aNormalSpeed = [&](double theU) {
    return RevolvedNormalSpeed(aSurface, theVelocity, theU, theV);
  };

  // Along a parallel the rotated normal and position yield only the
  // constant and first harmonic terms of the normal velocity.
  const double aZero       = aNormalSpeed(aFirst);
  const double aHalf       = aNormalSpeed(aFirst + aPi);
  const double aConstant   = 0.5 * (aZero + aHalf);
  const double aCosine     = 0.5 * (aZero - aHalf);
  const double aSine       = aNormalSpeed(aFirst + 0.5 * aPi) - aConstant;
  const double anAmplitude = std::hypot(aCosine, aSine);
  if (anAmplitude <= Precision::Confusion() || std::abs(aConstant) > anAmplitude)
    return {};

  const double aPhase     = aFirst + std::atan2(aSine, aCosine);
  const double aHalfWidth = std::acos(std::clamp(-aConstant / anAmplitude, -1.0, 1.0));
  std::array<std::optional<double>, 2> aResult;
  for (int i = 0; i < 2; ++i)
  {
    const double aBranch = aPhase + (i == 0 ? -aHalfWidth : aHalfWidth);
    const double aRoot =
      aBranch + aPeriod * std::ceil((aFirst - aBranch - Precision::Angular()) / aPeriod);
    if (aRoot > aLast + Precision::Angular())
      continue;

    const double                   aU = std::clamp(aRoot, aFirst, aLast);
    const BRepClass_FaceClassifier aClassifier(theFace,
                                               gp_Pnt2d(aU, theV),
                                               Precision::PConfusion());
    if (aClassifier.State() == TopAbs_IN || aClassifier.State() == TopAbs_ON)
      aResult[i] = aU;
  }
  return aResult;
}

inline std::array<std::optional<double>, 2> RevolvedScrewContactRoots(const TopoDS_Face& theFace,
                                                                      const gp_Ax1&      theAxis,
                                                                      double             theRise,
                                                                      double             theV)
{
  return RevolvedContactRoots(theFace, ContactVelocity::Screw(theAxis, theRise), theV);
}

inline std::vector<double> RevolvedScrewContactAngles(const TopoDS_Face& theFace,
                                                      const gp_Ax1&      theAxis,
                                                      double             theRise,
                                                      double             theV)
{
  std::vector<double> aResult;
  for (const auto& aRoot : RevolvedScrewContactRoots(theFace, theAxis, theRise, theV))
    if (aRoot)
      aResult.push_back(*aRoot);
  std::sort(aResult.begin(), aResult.end());
  if (aResult.size() == 2 && aResult[1] - aResult[0] <= Precision::Angular())
    aResult.pop_back();
  return aResult;
}

//! Locate normal-velocity crossings on the trimming wires. Sampling the
//! meridian alone can miss a whole contact branch in a narrow face sector.
inline std::vector<double> RevolvedTrimContactStations(const TopoDS_Face&         theFace,
                                                       const BRepAdaptor_Surface& theSurface,
                                                       const ContactVelocity&     theVelocity)
{
  std::vector<double> aStations{theSurface.FirstVParameter(), theSurface.LastVParameter()};
  for (TopExp_Explorer anIt(theFace, TopAbs_EDGE); anIt.More(); anIt.Next())
  {
    double     aFirst, aLast;
    const auto aPCurve =
      BRep_Tool::CurveOnSurface(TopoDS::Edge(anIt.Current()), theFace, aFirst, aLast);
    if (aPCurve.IsNull())
      continue;
    const auto speedAt = [&](double theParameter) {
      const auto aUV = aPCurve->Value(theParameter);
      gp_Pnt     aPoint;
      gp_Vec     aDu, aDv;
      theSurface.D1(aUV.X(), aUV.Y(), aPoint, aDu, aDv);
      return aDu.Crossed(aDv).Dot(theVelocity.At(aPoint));
    };
    constexpr int aSamples  = 32;
    double        aPrevious = speedAt(aFirst);
    for (int i = 1; i <= aSamples; ++i)
    {
      const double t        = aFirst + (aLast - aFirst) * i / aSamples;
      const double aCurrent = speedAt(t);
      if ((aCurrent < 0.0) != (aPrevious < 0.0))
      {
        double lo = aFirst + (aLast - aFirst) * (i - 1) / aSamples, hi = t;
        for (int j = 0; j < 48 && hi - lo > Precision::PConfusion(); ++j)
        {
          const double mid = 0.5 * (lo + hi);
          if ((speedAt(mid) < 0.0) == (aPrevious < 0.0))
            lo = mid;
          else
            hi = mid;
        }
        aStations.push_back(std::clamp(aPCurve->Value(0.5 * (lo + hi)).Y(),
                                       theSurface.FirstVParameter(),
                                       theSurface.LastVParameter()));
      }
      aPrevious = aCurrent;
    }
  }
  std::sort(aStations.begin(), aStations.end());
  return aStations;
}

//! Fit contact generators on trimmed analytic surfaces of revolution.
//! Contact birth and exit events are bracketed in the meridian parameter.
//! Surface fitting is checked between knots against the analytic roots.
inline std::vector<occ::handle<Geom_BSplineCurve>> FitRevolvedGenerators(
  const TopoDS_Face&     theFace,
  const ContactVelocity& theVelocity,
  double                 theTolerance)
{
  if (!std::isfinite(theTolerance) || theTolerance <= 0.0)
    throw Standard_ConstructionError("Invalid revolved contact fit tolerance");

  const BRepAdaptor_Surface aSurface(theFace);
  const double              aFirst = aSurface.FirstVParameter(), aLast = aSurface.LastVParameter();
  const double              aPi     = std::acos(-1.0);
  auto                      aProbes = RevolvedTrimContactStations(theFace, aSurface, theVelocity);
  const size_t              aTrimCount = aProbes.size();
  for (size_t i = 1; i < aTrimCount; ++i)
    aProbes.push_back(0.5 * (aProbes[i - 1] + aProbes[i]));
  constexpr int aScanCount = 128;
  for (int i = 1; i < aScanCount; ++i)
    aProbes.push_back(aFirst + (aLast - aFirst) * i / aScanCount);
  std::sort(aProbes.begin(), aProbes.end());
  aProbes.erase(std::unique(aProbes.begin(), aProbes.end()), aProbes.end());
  std::vector<occ::handle<Geom_BSplineCurve>> aGenerators;
  for (int branch = 0; branch < 2; ++branch)
  {
    const auto root = [&](double v) {
      return RevolvedContactRoots(theFace, theVelocity, v)[branch];
    };
    std::vector<double> aLimits{aFirst};
    bool                aPrevious = root(aFirst).has_value();
    for (size_t i = 1; i < aProbes.size(); ++i)
    {
      const double v        = aProbes[i];
      const bool   aPresent = root(v).has_value();
      if (aPresent != aPrevious)
      {
        double lo = aProbes[i - 1], hi = v;
        // A merging pair approaches its endpoint as a square root of the
        // meridian distance. Resolve the event down to machine precision.
        for (int j = 0; j < 64; ++j)
        {
          const double mid = (lo + hi) * 0.5;
          if (mid == lo || mid == hi)
            break;
          if (root(mid).has_value() == aPrevious)
            lo = mid;
          else
            hi = mid;
        }
        // Use the side where the branch is defined.
        aLimits.push_back(aPresent ? hi : lo);
      }
      aPrevious = aPresent;
    }
    aLimits.push_back(aLast);

    for (std::size_t i = 1; i < aLimits.size(); ++i)
    {
      const double a = aLimits[i - 1], b = aLimits[i];
      if (b - a <= Precision::PConfusion() || !root(0.5 * (a + b)))
        continue;
      const auto aStartU = root(a), anEndU = root(b);
      if (!aStartU || !anEndU)
        throw Standard_ConstructionError("Could not evaluate revolved contact ends");
      const double aUSpan     = std::abs(*anEndU - *aStartU);
      const bool   useAzimuth = b - a < 0.01 * aUSpan && aUSpan < aPi;
      const auto   aPointAt   = [&](double u) {
        if (useAzimuth)
        {
          if (u <= 0.0)
            return aSurface.Value(*aStartU, a);
          if (u >= 1.0)
            return aSurface.Value(*anEndU, b);
          const double aU = *aStartU + u * (*anEndU - *aStartU);
          double       lo = a, hi = b;
          const double aSpeedLo = RevolvedNormalSpeed(aSurface, theVelocity, aU, lo);
          const double aSpeedHi = RevolvedNormalSpeed(aSurface, theVelocity, aU, hi);
          if ((aSpeedLo < 0.0) == (aSpeedHi < 0.0))
            throw Standard_ConstructionError("Could not bracket near-parallel revolved contact");
          for (int i = 0; i < 48 && hi - lo > Precision::PConfusion(); ++i)
          {
            const double mid = 0.5 * (lo + hi);
            if ((RevolvedNormalSpeed(aSurface, theVelocity, aU, mid) < 0.0) == (aSpeedLo < 0.0))
              lo = mid;
            else
              hi = mid;
          }
          return aSurface.Value(aU, 0.5 * (lo + hi));
        }
        const double v  = u <= 0.0   ? a
                          : u >= 1.0 ? b
                                     : a + (b - a) * (0.5 * (1.0 - std::cos(aPi * u)));
        const auto   aU = root(v);
        if (!aU)
          throw Standard_ConstructionError("Revolved contact branch left supporting face");
        return aSurface.Value(*aU, v);
      };

      occ::handle<Geom_BSplineCurve> aCurve;
      for (int count = 17; count <= 257; count = 2 * count - 1)
      {
        const auto aPoints =
          occ::handle<NCollection_HArray1<gp_Pnt>>(new NCollection_HArray1<gp_Pnt>(1, count));
        const auto aParameters =
          occ::handle<NCollection_HArray1<double>>(new NCollection_HArray1<double>(1, count));
        for (int j = 1; j <= count; ++j)
        {
          const double u = double(j - 1) / (count - 1);
          aParameters->SetValue(j, u);
          aPoints->SetValue(j, aPointAt(u));
        }
        GeomAPI_Interpolate aFit(aPoints, aParameters, false, Precision::Confusion());
        aFit.Perform();
        if (!aFit.IsDone())
          continue;

        double anError = 0.0;
        for (int j = 1; j < 4 * (count - 1); ++j)
        {
          const double u = double(j) / (4 * (count - 1));
          anError        = std::max(anError, aFit.Curve()->Value(u).Distance(aPointAt(u)));
        }
        if (anError <= theTolerance)
        {
          aCurve = aFit.Curve();
          break;
        }
      }
      if (aCurve.IsNull())
        throw Standard_ConstructionError("Revolved contact generator missed the fit tolerance");
      aGenerators.push_back(aCurve);
    }
  }
  return aGenerators;
}

inline std::vector<occ::handle<Geom_BSplineCurve>> FitRevolvedScrewGenerators(
  const TopoDS_Face& theFace,
  const gp_Ax1&      theAxis,
  double             theRise,
  double             theTolerance)
{
  return FitRevolvedGenerators(theFace, ContactVelocity::Screw(theAxis, theRise), theTolerance);
}

//! Candidate spans on a manifold sharp edge where the two incident face
//! normals have opposite velocities. Parameters belong to the 3D edge.
//! The caller must still reject locally concave and globally buried spans.
inline std::vector<std::pair<double, double>> SharpContactSpans(const TopoDS_Edge& theEdge,
                                                                const TopoDS_Face& theFirstFace,
                                                                const TopoDS_Face& theSecondFace,
                                                                const ContactVelocity& theVelocity)
{
  if (theFirstFace.IsSame(theSecondFace) || BRep_Tool::Degenerated(theEdge))
    return {};
  if (!BRep_Tool::SameParameter(theEdge))
    throw Standard_ConstructionError("Screw contact edge requires same-parameter pcurves");

  const BRepAdaptor_Curve   anEdge(theEdge);
  const BRepAdaptor_Surface aFirstSurface(theFirstFace), aSecondSurface(theSecondFace);
  double                    aFirst, aLast;
  const auto aFirstPCurve  = BRep_Tool::CurveOnSurface(theEdge, theFirstFace, aFirst, aLast);
  const auto aSecondPCurve = BRep_Tool::CurveOnSurface(theEdge, theSecondFace, aFirst, aLast);
  if (aFirstPCurve.IsNull() || aSecondPCurve.IsNull())
    throw Standard_ConstructionError("Missing screw contact pcurve");

  aFirst                  = anEdge.FirstParameter();
  aLast                   = anEdge.LastParameter();
  const auto aNormalSpeed = [&](bool theOnFirst, double theT) {
    const auto&    aFace    = theOnFirst ? theFirstFace : theSecondFace;
    const auto&    aSurface = theOnFirst ? aFirstSurface : aSecondSurface;
    const auto&    aPCurve  = theOnFirst ? aFirstPCurve : aSecondPCurve;
    const gp_Pnt2d aUV      = aPCurve->Value(theT);
    gp_Pnt         aPoint;
    gp_Vec         aDu, aDv;
    aSurface.D1(aUV.X(), aUV.Y(), aPoint, aDu, aDv);
    gp_Vec aNormal = aDu.Crossed(aDv);
    if (aNormal.SquareMagnitude() <= gp::Resolution())
      return 0.0;
    aNormal.Normalize();
    if (aFace.Orientation() == TopAbs_REVERSED)
      aNormal.Reverse();
    return aNormal.Dot(theVelocity.At(anEdge.Value(theT)));
  };

  std::vector<double> aBreaks{aFirst, aLast};
  constexpr int       aSamples = 64;
  for (const bool anOnFirst : {true, false})
  {
    double aPrevious = aNormalSpeed(anOnFirst, aFirst);
    for (int i = 1; i <= aSamples; ++i)
    {
      const double u        = aFirst + (aLast - aFirst) * i / aSamples;
      const double aCurrent = aNormalSpeed(anOnFirst, u);
      if ((aPrevious < 0.0) != (aCurrent < 0.0))
      {
        double lo = aFirst + (aLast - aFirst) * (i - 1) / aSamples, hi = u;
        for (int j = 0; j < 48 && hi - lo > Precision::PConfusion(); ++j)
        {
          const double mid = (lo + hi) * 0.5;
          if ((aNormalSpeed(anOnFirst, mid) < 0.0) == (aPrevious < 0.0))
            lo = mid;
          else
            hi = mid;
        }
        aBreaks.push_back((lo + hi) * 0.5);
      }
      aPrevious = aCurrent;
    }
  }

  std::sort(aBreaks.begin(), aBreaks.end());
  std::vector<std::pair<double, double>> aResult;
  for (std::size_t i = 1; i < aBreaks.size(); ++i)
  {
    const double lo = aBreaks[i - 1], hi = aBreaks[i];
    if (hi - lo <= Precision::PConfusion())
      continue;
    const double mid = (lo + hi) * 0.5;
    const double a = aNormalSpeed(true, mid), b = aNormalSpeed(false, mid);
    if (std::min(a, b) >= -Precision::Confusion() || std::max(a, b) <= Precision::Confusion())
      continue;

    if (!aResult.empty() && lo - aResult.back().second <= Precision::PConfusion())
      aResult.back().second = hi;
    else
      aResult.emplace_back(lo, hi);
  }
  return aResult;
}

inline std::vector<std::pair<double, double>> SharpScrewContactSpans(
  const TopoDS_Edge& theEdge,
  const TopoDS_Face& theFirstFace,
  const TopoDS_Face& theSecondFace,
  const gp_Ax1&      theAxis,
  double             theRise)
{
  return SharpContactSpans(theEdge,
                           theFirstFace,
                           theSecondFace,
                           ContactVelocity::Screw(theAxis, theRise));
}

//! Sweep a contact generator continuously along a screw motion. The
//! generator's rational control polygon and weights are retained exactly.
//! A cubic Hermite fit in the motion parameter bounds the approximation by
//! radius * angularStep^4 / 384, using the largest control-point radius.
inline occ::handle<Geom_BSplineSurface> FitScrewContactRibbon(
  const occ::handle<Geom_Curve>& theGenerator,
  const gp_Ax1&                  theAxis,
  double                         theRise,
  double                         theAngle,
  double                         theTolerance)
{
  if (theGenerator.IsNull() || !std::isfinite(theRise) || !std::isfinite(theAngle)
      || theAngle <= Precision::Angular() || !std::isfinite(theTolerance) || theTolerance <= 0.0)
    throw Standard_ConstructionError("Invalid screw contact ribbon input");

  const auto aGenerator = GeomConvert::CurveToBSplineCurve(theGenerator);
  if (aGenerator->IsPeriodic())
    aGenerator->SetNotPeriodic();
  const gp_Vec anOmega(theAxis.Direction());
  double       aRadius = 0.0;
  for (int i = 1; i <= aGenerator->NbPoles(); ++i)
  {
    const gp_Vec aDisplacement(theAxis.Location(), aGenerator->Pole(i));
    aRadius = std::max(aRadius, anOmega.Crossed(aDisplacement).Magnitude());
  }

  const double aPi = std::acos(-1.0);
  const double aStepLimit =
    aRadius == 0.0 ? aPi / 4.0
                   : std::min(aPi / 4.0, std::sqrt(std::sqrt((384.0 * theTolerance) / aRadius)));
  const double aRequired = std::ceil(theAngle / aStepLimit);
  if (!(aRequired <= 4096.0))
    throw Standard_ConstructionError("Screw contact ribbon exceeded fitting capacity");
  const int    aCount = std::max(1, int(aRequired));
  const double aStep  = theAngle / aCount;

  NCollection_Array1<double> aKnots(1, aCount + 1);
  NCollection_Array1<int>    aMults(1, aCount + 1);
  for (int j = 0; j <= aCount; ++j)
  {
    aKnots(j + 1) = j == aCount ? theAngle : j * aStep;
    aMults(j + 1) = j == 0 || j == aCount ? 4 : 3;
  }

  NCollection_Array2<gp_Pnt> aPoles(1, aGenerator->NbPoles(), 1, 3 * aCount + 1);
  NCollection_Array2<double> aWeights(1, aGenerator->NbPoles(), 1, 3 * aCount + 1);
  for (int i = 1; i <= aGenerator->NbPoles(); ++i)
  {
    const gp_Pnt aPole     = aGenerator->Pole(i);
    gp_Pnt       aPrevious = aPole;
    const auto   aVelocity = [&](const gp_Pnt& thePoint) {
      return anOmega.Crossed(gp_Vec(theAxis.Location(), thePoint)) + theRise * anOmega;
    };
    aPoles(i, 1)   = aPole;
    aWeights(i, 1) = aGenerator->Weight(i);
    for (int j = 1; j <= aCount; ++j)
    {
      const double t         = aKnots(j + 1);
      const gp_Pnt aPosition = aPole.Rotated(theAxis, t).Translated((theRise * t) * anOmega);
      aPoles(i, 3 * j - 1)   = aPrevious.Translated((aStep / 3.0) * aVelocity(aPrevious));
      aPoles(i, 3 * j)       = aPosition.Translated((-aStep / 3.0) * aVelocity(aPosition));
      aPoles(i, 3 * j + 1)   = aPosition;
      for (int k = 3 * j - 1; k <= 3 * j + 1; ++k)
        aWeights(i, k) = aGenerator->Weight(i);
      aPrevious = aPosition;
    }
  }

  return new Geom_BSplineSurface(aPoles,
                                 aWeights,
                                 aGenerator->Knots(),
                                 aKnots,
                                 aGenerator->Multiplicities(),
                                 aMults,
                                 aGenerator->Degree(),
                                 3);
}

//! Retain the part of a planar contact chord without a strict local fold.
//! Both endpoints must satisfy the contact condition n dot velocity = 0.
//! The signed plane distance along the inverse orbit has second derivative
//! n dot (w cross (w cross (p - origin))). Axial advance has zero second
//! derivative. Negative distance is inside the solid, for an outward normal.
//! Zero is kept for higher-order and finite-motion classification.
inline std::optional<OrbitInterval> ExposedPlanarChord(const gp_Pnt& theFirst,
                                                       const gp_Pnt& theLast,
                                                       const gp_Dir& theOutward,
                                                       const gp_Ax1& theAxis)
{
  const gp_Vec anOmega(theAxis.Direction());
  const gp_Vec aNormal(theOutward);
  const auto   aSecondDerivative = [&](const gp_Pnt& thePoint) {
    return aNormal.Dot(anOmega.Crossed(anOmega.Crossed(gp_Vec(theAxis.Location(), thePoint))));
  };
  const double aFirst = aSecondDerivative(theFirst);
  const double aLast  = aSecondDerivative(theLast);
  if (aFirst >= 0.0 && aLast >= 0.0)
    return OrbitInterval{0.0, 1.0};
  if (aFirst <= 0.0 && aLast <= 0.0)
    return std::nullopt;
  const double aRoot = aFirst / (aFirst - aLast);
  return aFirst > 0.0 ? OrbitInterval{0.0, aRoot} : OrbitInterval{aRoot, 1.0};
}

//! Boundary of the swept occupancy and the stationary interval that supplied it.
struct SweptOrbitInterval
{
  double      First;
  double      Last;
  std::size_t FirstSource;
  std::size_t LastSource;
};

//! Dilate stationary occupancy by [0, theAngle] and unite overlapping intervals.
//! Input intervals are finite and ordered internally (First <= Last).
//! A positive angle advances the upper bound. A negative angle advances the lower bound.
//! Angular coordinates must be unwrapped, including on orbits spanning several turns.
inline std::vector<SweptOrbitInterval> SweepOrbitIntervals(
  const std::vector<OrbitInterval>& theOccupancy,
  double                            theAngle)
{
  const double aFirstShift = std::min(0.0, theAngle);
  const double aLastShift  = std::max(0.0, theAngle);

  std::vector<SweptOrbitInterval> aResult;
  aResult.reserve(theOccupancy.size());
  for (std::size_t anIndex = 0; anIndex < theOccupancy.size(); ++anIndex)
  {
    const OrbitInterval& anInterval = theOccupancy[anIndex];
    aResult.push_back(
      {anInterval.First + aFirstShift, anInterval.Last + aLastShift, anIndex, anIndex});
  }
  std::stable_sort(aResult.begin(),
                   aResult.end(),
                   [](const SweptOrbitInterval& theLeft, const SweptOrbitInterval& theRight) {
                     return theLeft.First < theRight.First;
                   });

  std::size_t aCount = 0;
  for (const SweptOrbitInterval& anInterval : aResult)
  {
    if (aCount == 0 || anInterval.First > aResult[aCount - 1].Last)
    {
      aResult[aCount++] = anInterval;
      continue;
    }

    SweptOrbitInterval& aMerged = aResult[aCount - 1];
    if (anInterval.Last > aMerged.Last)
    {
      aMerged.Last       = anInterval.Last;
      aMerged.LastSource = anInterval.LastSource;
    }
  }
  aResult.resize(aCount);
  return aResult;
}

//! Coordinates for nonzero-rise screw motion. The planar key stays constant
//! along each orbit. The axial coordinate lifts a point on that orbit back
//! into three-dimensional space without wrapping turns.
class ScrewOrbitCoordinates
{
public:
  ScrewOrbitCoordinates(const gp_Ax1& theAxis, double theRise)
      : myFrame(theAxis.Location(), theAxis.Direction()),
        myRise(theRise)
  {
    if (!std::isfinite(myRise) || myRise == 0.0)
      throw Standard_ConstructionError("Screw orbit coordinates require nonzero rise");
  }

  double Axial(const gp_Pnt& thePoint) const
  {
    return gp_Vec(myFrame.Location(), thePoint).Dot(gp_Vec(myFrame.Direction()));
  }

  gp_Pnt2d Key(const gp_Pnt& thePoint) const
  {
    const double anAxial = Axial(thePoint);
    const gp_Vec aDisplacement(myFrame.Location(), thePoint);
    const double aPhase = Phase(anAxial);
    const double aCos = std::cos(aPhase), aSin = std::sin(aPhase);
    const double anX = aDisplacement.Dot(gp_Vec(myFrame.XDirection()));
    const double anY = aDisplacement.Dot(gp_Vec(myFrame.YDirection()));
    return gp_Pnt2d(aCos * anX + aSin * anY, -aSin * anX + aCos * anY);
  }

  gp_Pnt Lift(const gp_Pnt2d& theKey, double theAxial) const
  {
    const double aPhase = Phase(theAxial);
    const double aCos = std::cos(aPhase), aSin = std::sin(aPhase);
    return myFrame.Location().Translated(
      (aCos * theKey.X() - aSin * theKey.Y()) * gp_Vec(myFrame.XDirection())
      + (aSin * theKey.X() + aCos * theKey.Y()) * gp_Vec(myFrame.YDirection())
      + theAxial * gp_Vec(myFrame.Direction()));
  }

  //! Angular separation of points known to have the same planar key.
  double Lag(const gp_Pnt& theFrom, const gp_Pnt& theTo) const
  {
    return gp_Vec(theFrom, theTo).Dot(gp_Vec(myFrame.Direction())) / myRise;
  }

private:
  double Phase(double theAxial) const
  {
    return std::remainder(theAxial / myRise, 2.0 * std::acos(-1.0));
  }

  gp_Ax2 myFrame;
  double myRise;
};

//! Intersect the stationary solid with an exact screw orbit in an angular window.
class OrbitIntersector
{
public:
  OrbitIntersector(const TopoDS_Solid& theSolid, double theTolerance)
      : myClassifier(theSolid),
        myTolerance(theTolerance)
  {
    myIntersector.Load(theSolid, theTolerance);
  }

  std::vector<OrbitInterval> Occupied(const gp_Pnt& thePoint,
                                      const gp_Ax1& theAxis,
                                      double        theRise,
                                      double        theFirst,
                                      double        theLast)
  {
    if (!(theFirst < theLast) || !std::isfinite(theFirst) || !std::isfinite(theLast)
        || !std::isfinite(theRise))
    {
      throw Standard_ConstructionError("Invalid screw orbit range");
    }

    const gp_Vec aDirection(theAxis.Direction());
    const gp_Vec anOffset(theAxis.Location(), thePoint);
    const gp_Vec aRadial = anOffset - anOffset.Dot(aDirection) * aDirection;
    const double aRadius = aRadial.Magnitude();
    if (aRadius <= Precision::Confusion() && std::abs(theRise) <= Precision::Confusion())
    {
      myClassifier.Perform(thePoint, myTolerance);
      return myClassifier.State() == TopAbs_OUT ? std::vector<OrbitInterval>{}
                                                : std::vector<OrbitInterval>{{theFirst, theLast}};
    }

    // On a coaxial cylinder the screw orbit is a straight line in (angle, height).
    occ::handle<Adaptor3d_Curve> aCurve;
    double                       aScale;
    if (aRadius <= Precision::Confusion())
    {
      aScale = std::abs(theRise);
      const auto aLine =
        occ::handle<Geom_Line>(new Geom_Line(thePoint, gp_Dir(theRise * aDirection)));
      aCurve = new GeomAdaptor_Curve(aLine, theFirst * aScale, theLast * aScale);
    }
    else
    {
      aScale                 = std::hypot(1.0, theRise);
      const gp_Pnt anOrigin  = thePoint.Translated(-aRadial);
      const auto   aCylinder = occ::handle<Geom_CylindricalSurface>(
        new Geom_CylindricalSurface(gp_Ax3(anOrigin, theAxis.Direction(), gp_Dir(aRadial)),
                                    aRadius));
      const auto aLine =
        occ::handle<Geom2d_Line>(new Geom2d_Line(gp_Pnt2d(), gp_Dir2d(1.0, theRise)));
      aCurve = new Adaptor3d_CurveOnSurface(
        new Geom2dAdaptor_Curve(aLine, theFirst * aScale, theLast * aScale),
        new GeomAdaptor_Surface(aCylinder));
    }

    std::vector<double>        aBreaks{theFirst, theLast};
    std::vector<OrbitInterval> anOccupied;
    const double               aParamTolerance = myTolerance / std::hypot(aRadius, theRise);
    const double               aQuarterTurn    = std::acos(-1.0) / 2.0;
    for (double aStart = theFirst; aStart < theLast;)
    {
      const double anEnd = std::min(aStart + aQuarterTurn, theLast);
      if (!(anEnd > aStart))
        throw Standard_ConstructionError("Screw orbit range exceeds parameter precision");
      const auto aSegment = aCurve->Trim(aStart * aScale, anEnd * aScale, Precision::PConfusion());
      myIntersector.Perform(aSegment, aStart * aScale, anEnd * aScale);
      if (!myIntersector.IsDone())
        throw Standard_ConstructionError("Could not intersect screw orbit with stationary solid");
      for (int anIndex = 1; anIndex <= myIntersector.NbPnt(); ++anIndex)
      {
        const double aHit = std::clamp(myIntersector.WParameter(anIndex) / aScale, aStart, anEnd);
        aBreaks.push_back(aHit);
        anOccupied.push_back({aHit, aHit});
      }
      aStart = anEnd;
    }

    std::sort(aBreaks.begin(), aBreaks.end());
    for (std::size_t anIndex = 1; anIndex < aBreaks.size(); ++anIndex)
    {
      const double a = aBreaks[anIndex - 1], b = aBreaks[anIndex];
      if (b - a <= aParamTolerance)
        continue;
      myClassifier.Perform(aCurve->Value(aScale * (a + (b - a) / 2.0)), myTolerance);
      if (myClassifier.State() == TopAbs_IN || myClassifier.State() == TopAbs_ON)
        anOccupied.push_back({a, b});
    }

    const auto aMerged = SweepOrbitIntervals(anOccupied, 0.0);
    anOccupied.clear();
    for (const auto& anInterval : aMerged)
      anOccupied.push_back({anInterval.First, anInterval.Last});
    return anOccupied;
  }

private:
  IntCurvesFace_ShapeIntersector myIntersector;
  BRepClass3d_SolidClassifier    myClassifier;
  double                         myTolerance;
};
} // namespace BRepSolidSweep

#endif

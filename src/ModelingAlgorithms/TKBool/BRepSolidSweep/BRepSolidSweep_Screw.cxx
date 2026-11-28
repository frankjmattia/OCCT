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

#include "BRepSolidSweep_ChordFit.hxx"
#include "BRepSolidSweep_OrbitIntervals.hxx"
#include <BRepSolidSweep_Sweep.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <Precision.hxx>
#include <Standard_ConstructionError.hxx>
#include <TopExp_Explorer.hxx>
#include <gp_Cone.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Pln.hxx>
#include <gp_Torus.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepLib.hxx>
#include <GC_MakeSegment2d.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <Geom_Curve.hxx>
#include <Geom_SurfaceOfRevolution.hxx>
#include <NCollection_HArray1.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRep_Tool.hxx>
#include <Geom_BSplineSurface.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <Standard_Failure.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Ax1.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Builder.hxx>
#include <GeomAPI.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomProjLib.hxx>
#include <Geom_Plane.hxx>
#include <GProp_GProps.hxx>
#include <NCollection_HSequence.hxx>
#include <ShapeAnalysis_FreeBounds.hxx>
#include <TopExp.hxx>
#include <string_view>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace
{
BRepSolidSweep::ScrewStationSchedule gScrewSchedule =
  BRepSolidSweep::ScrewStationSchedule::Curvature;
std::vector<BRepSolidSweep::ScrewPatchFit> gScrewFits;
} // namespace

void BRepSolidSweep::SetScrewStationSchedule(ScrewStationSchedule theSchedule)
{
  gScrewSchedule = theSchedule;
}

BRepSolidSweep::ScrewStationSchedule BRepSolidSweep::ActiveScrewStationSchedule()
{
  return gScrewSchedule;
}

void BRepSolidSweep::ClearScrewPatchFits()
{
  gScrewFits.clear();
}

const std::vector<BRepSolidSweep::ScrewPatchFit>& BRepSolidSweep::ScrewPatchFits()
{
  return gScrewFits;
}

namespace
{
const double pi = std::acos(-1.0);

struct Band
{
  TopoDS_Face Face;
  gp_Ax1      Axis;
  double      Radius;
};

struct AxisymmetricBody
{
  gp_Pnt            Base;
  gp_Dir            Axis;
  double            Width;
  double            BackRadius;
  double            FrontRadius;
  double            MaxRadius;
  std::vector<Band> Bands;
};

gp_Circ readDisk(const TopoDS_Face& theFace)
{
  TopExp_Explorer edges(theFace, TopAbs_EDGE);
  if (!edges.More())
  {
    throw Standard_ConstructionError("Expected a circular end cap");
  }

  const BRepAdaptor_Curve curve(TopoDS::Edge(edges.Current()));

  edges.Next();
  if (edges.More() || curve.GetType() != GeomAbs_Circle
      || std::abs(curve.LastParameter() - curve.FirstParameter() - 2.0 * pi) > Precision::Angular())
  {
    throw Standard_ConstructionError("Expected a disk cap without a bore");
  }

  return curve.Circle();
}

AxisymmetricBody readBody(const TopoDS_Solid&                      theSolid,
                          const BRepSolidSweep_Sweep::ScrewMotion& theMotion)
{
  if (theSolid.IsNull() || !BRepCheck_Analyzer(theSolid).IsValid())
  {
    throw Standard_ConstructionError("Expected a valid axisymmetric solid");
  }

  std::vector<Band>    bands;
  std::vector<gp_Circ> caps;
  for (TopExp_Explorer faces(theSolid, TopAbs_FACE); faces.More(); faces.Next())
  {
    const TopoDS_Face&        aFace = TopoDS::Face(faces.Current());
    const BRepAdaptor_Surface surface(aFace);
    if (surface.GetType() == GeomAbs_Plane)
    {
      caps.push_back(readDisk(aFace));
      continue;
    }

    gp_Ax1 axis;
    double radius;
    switch (surface.GetType())
    {
      case GeomAbs_Cylinder:
        axis   = surface.Cylinder().Axis();
        radius = surface.Cylinder().Radius();
        break;
      case GeomAbs_Cone:
        axis   = surface.Cone().Axis();
        radius = surface.Cone().RefRadius();
        break;
      case GeomAbs_Torus:
        axis   = surface.Torus().Axis();
        radius = surface.Torus().MajorRadius() + surface.Torus().MinorRadius();
        break;
      default:
        throw Standard_ConstructionError("Unsupported side face on axisymmetric solid");
    }

    if (std::abs(surface.LastUParameter() - surface.FirstUParameter() - 2.0 * pi)
        > Precision::Angular())
    {
      throw Standard_ConstructionError("Expected complete revolved side faces");
    }

    bands.push_back({aFace, axis, radius});
  }
  if (bands.empty() || caps.size() != 2)
  {
    throw Standard_ConstructionError("Expected revolved side faces and two disk caps");
  }

  gp_Vec along(caps.front().Location(), caps.back().Location());
  if (along.Magnitude() <= Precision::Confusion())
  {
    throw Standard_ConstructionError("End caps must have distinct centers");
  }
  const gp_Vec vertical(theMotion.Axis.Direction());
  const gp_Vec velocity =
    vertical.Crossed(gp_Vec(theMotion.Axis.Location(), caps.front().Location()))
    + (theMotion.AdvancePerTurn / (2.0 * pi)) * vertical;
  if (along.Dot(velocity) < 0.0)
  {
    std::swap(caps.front(), caps.back());
    along.Reverse();
  }

  const gp_Dir  axis(along);
  const gp_Pnt& base = caps.front().Location();
  const gp_Vec  direction(axis);
  for (const gp_Circ& cap : caps)
  {
    if (!cap.Axis().Direction().IsParallel(axis, Precision::Angular()))
    {
      throw Standard_ConstructionError("End caps must be perpendicular to the body axis");
    }
  }
  double maxRadius = std::max(caps.front().Radius(), caps.back().Radius());
  for (const Band& band : bands)
  {
    const gp_Vec offset(base, band.Axis.Location());
    if (!band.Axis.Direction().IsParallel(axis, Precision::Angular())
        || offset.Crossed(direction).Magnitude() > 10.0 * Precision::Confusion())
    {
      throw Standard_ConstructionError("Revolved side faces must be coaxial");
    }
    maxRadius = std::max(maxRadius, band.Radius);
  }
  return {base,
          axis,
          along.Magnitude(),
          caps.front().Radius(),
          caps.back().Radius(),
          maxRadius,
          std::move(bands)};
}

// In the co-rotating frame the velocity field is stationary. At each
// meridian station the tangency equation is a*cos(phi) + b*sin(phi) + c = 0.
struct ScrewEnvelope
{
  struct MeridianBand
  {
    occ::handle<BRepAdaptor_Surface> Surface;
    occ::handle<Geom_Curve>          Curve;
    double                           U;
    bool                             Reversed;
  };

  gp_Pnt                    base;
  gp_Vec                    axis;
  double                    width, backFaceRadius, frontFaceRadius, maxRadius;
  gp_Ax1                    rotation;
  gp_Vec                    vertical, radialX, radialY;
  double                    rise, angle;
  int                       timeIntervals;
  std::vector<MeridianBand> bands;

  ScrewEnvelope(const AxisymmetricBody& body, const BRepSolidSweep_Sweep::ScrewMotion& motion)
      : base(body.Base),
        axis(body.Axis),
        width(body.Width),
        backFaceRadius(body.BackRadius),
        frontFaceRadius(body.FrontRadius),
        maxRadius(body.MaxRadius),
        rotation(motion.Axis),
        vertical(motion.Axis.Direction()),
        rise(motion.AdvancePerTurn / (2.0 * pi)),
        angle(motion.Angle),
        timeIntervals(static_cast<int>(std::ceil(angle / (pi / 64.0))))
  {
    const gp_Vec offset(rotation.Location(), base);
    const gp_Vec baseRadial = offset - offset.Dot(vertical) * vertical;
    const gp_Vec projected  = baseRadial - baseRadial.Dot(axis) * axis;
    radialX                 = projected.Magnitude() > Precision::Confusion()
                                ? projected.Normalized()
                                : gp_Vec(gp_Ax2(base, body.Axis).XDirection());
    radialY                 = axis.Crossed(radialX);

    for (const Band& band : body.Bands)
    {
      const occ::handle<BRepAdaptor_Surface> surface = new BRepAdaptor_Surface(band.Face);
      const double                           u       = surface->FirstUParameter() + pi / 4.0;
      const double v = 0.5 * (surface->FirstVParameter() + surface->LastVParameter());
      const gp_Vec fromBase(base, surface->Value(u, v));
      const gp_Vec fromAxis                  = fromBase - fromBase.Dot(axis) * axis;
      const occ::handle<Geom_Curve> meridian = BRep_Tool::Surface(band.Face)->UIso(u);
      meridian->Rotate(gp_Ax1(base, body.Axis), fromAxis.AngleWithRef(radialX, axis));
      bands.push_back({surface, meridian, u, band.Face.Orientation() == TopAbs_REVERSED});
    }
  }

  struct Sample
  {
    gp_Pnt point;
    gp_Vec derivative;
  };

  struct Meridian
  {
    double z, r, nz, nr;
  };

  double meridianParameter(int part, double u) const
  {
    const auto& surface = bands[part].Surface;
    return surface->FirstVParameter()
           + u * (surface->LastVParameter() - surface->FirstVParameter());
  }

  Meridian profile(int part, double u) const
  {
    const MeridianBand& band = bands[part];
    gp_Pnt              point;
    gp_Vec              du, dv;
    band.Surface->D1(band.U, meridianParameter(part, u), point, du, dv);
    const gp_Vec normal = (band.Reversed ? -1.0 : 1.0) * du.Crossed(dv).Normalized();
    const gp_Vec fromBase(base, point);
    const double z      = fromBase.Dot(axis);
    const gp_Vec radial = fromBase - z * axis;
    return {z, radial.Magnitude(), normal.Dot(axis), normal.Dot(radial.Normalized())};
  }

  void validateMeridian()
  {
    struct Bounds
    {
      Meridian First;
      Meridian Last;
    };

    std::vector<Bounds> bounds;
    for (int part = 0; part < static_cast<int>(bands.size()); ++part)
    {
      Meridian   first    = profile(part, 0.0);
      Meridian   last     = profile(part, 1.0);
      const bool reversed = first.z > last.z;
      if (reversed)
      {
        std::swap(first, last);
      }
      double lastZ = first.z;
      for (int i = 0; i <= 32; ++i)
      {
        const double   u      = i / 32.0;
        const Meridian sample = profile(part, reversed ? 1.0 - u : u);
        if (sample.nr < -1.0e-7 || sample.z < lastZ - Precision::Confusion())
        {
          throw Standard_ConstructionError("Expected an outward, axially monotone meridian");
        }
        maxRadius = std::max(maxRadius, sample.r);
        lastZ     = sample.z;
      }
      bounds.push_back({first, last});
    }
    std::sort(bounds.begin(), bounds.end(), [](const Bounds& left, const Bounds& right) {
      return left.First.z < right.First.z;
    });
    Meridian   previous{0.0, backFaceRadius, -1.0, 0.0};
    const auto requireTangent = [](const Meridian& first, const Meridian& second) {
      if (std::hypot(first.z - second.z, first.r - second.r) > 10.0 * Precision::Confusion()
          || std::hypot(first.nz - second.nz, first.nr - second.nr) > 1.0e-6)
      {
        throw Standard_ConstructionError("Screw sweep requires tangent joins in the meridian");
      }
    };
    for (const Bounds& band : bounds)
    {
      requireTangent(previous, band.First);
      previous = band.Last;
    }
    requireTangent(previous, {width, frontFaceRadius, 1.0, 0.0});
  }

  struct ContactEquation
  {
    double a, b, c;
  };

  ContactEquation equation(const Meridian& p)
  {
    const gp_Vec reference = vertical.Crossed(gp_Vec(rotation.Location(), base)) + rise * vertical;
    const gp_Vec atStation = reference + p.z * vertical.Crossed(axis);
    return {p.nr * radialX.Dot(atStation) + p.nz * p.r * axis.Dot(vertical.Crossed(radialX)),
            p.nr * radialY.Dot(atStation) + p.nz * p.r * axis.Dot(vertical.Crossed(radialY)),
            p.nz * axis.Dot(reference)};
  }

  double discriminant(int part, double u)
  {
    const auto e = equation(profile(part, u));
    return std::hypot(e.a, e.b) - std::abs(e.c);
  }

  double root(int part, double a, double b)
  {
    const bool startPositive = discriminant(part, a) >= 0.0;
    for (int i = 0; i < 50; ++i)
    {
      const double mid = (a + b) * 0.5;
      if ((discriminant(part, mid) >= 0.0) == startPositive)
        a = mid;
      else
        b = mid;
    }
    return (a + b) * 0.5;
  }

  gp_Pnt contact(int part, double a, double b, int branch, double u)
  {
    // Ease into event stations so the double root has a regular parameter.
    const double   eased = 0.5 * (1.0 - std::cos(pi * u));
    const Meridian p     = profile(part, a + (b - a) * eased);
    const auto     e     = equation(p);
    const double   phase = std::atan2(e.b, e.a);
    const double   range = std::acos(std::clamp(-e.c / std::hypot(e.a, e.b), -1.0, 1.0));
    const double   phi   = phase + branch * range;
    return base.Translated(p.z * axis + p.r * (std::cos(phi) * radialX + std::sin(phi) * radialY));
  }

  struct EndpointClearance
  {
    double distance;
    int    support;
  };

  EndpointClearance endpointClearance(const gp_Pnt& point) const
  {
    const gp_Vec      local(base, point);
    const double      z = local.Dot(axis);
    const double      r = (local - z * axis).Magnitude();
    EndpointClearance result =
      z < width - z ? EndpointClearance{-z, -1} : EndpointClearance{z - width, -2};
    const double section = std::clamp(z, 0.0, width);
    for (int part = 0; part < int(bands.size()); ++part)
    {
      const double start = profile(part, 0.0).z, end = profile(part, 1.0).z;
      if (section < std::min(start, end) - Precision::Confusion()
          || section > std::max(start, end) + Precision::Confusion())
        continue;
      const bool reversed = start > end;
      double     lo = 0.0, hi = 1.0;
      for (int i = 0; i < 48; ++i)
      {
        const double mid = (lo + hi) / 2.0;
        if (profile(part, reversed ? 1.0 - mid : mid).z < section)
          lo = mid;
        else
          hi = mid;
      }
      const double u      = (lo + hi) / 2.0;
      const double radial = r - profile(part, reversed ? 1.0 - u : u).r;
      if (radial > result.distance)
        result = {radial, part};
      return result;
    }
    throw Standard_Failure("Stationary section missing from endpoint meridian");
  }

  double localContactClearance(int part, double first, double last, int branch, double u)
  {
    const auto&  band    = bands[part];
    const double station = first + (last - first) * (0.5 * (1.0 - std::cos(pi * u)));
    const auto   point   = contact(part, first, last, branch, u);
    gp_Pnt       reference;
    gp_Vec       du, dv, duu, dvv, duv;
    band.Surface->D2(band.U, meridianParameter(part, station), reference, du, dv, duu, dvv, duv);
    const gp_Vec fromBase(base, reference);
    const gp_Pnt onAxis  = base.Translated(fromBase.Dot(axis) * axis);
    const double azimuth = gp_Vec(onAxis, reference).AngleWithRef(gp_Vec(onAxis, point), axis);
    const gp_Ax1 bodyAxis(base, gp_Dir(axis));
    for (auto* derivative : {&du, &dv, &duu, &duv, &dvv})
      derivative->Rotate(bodyAxis, azimuth);

    const gp_Dir outward((band.Reversed ? -1.0 : 1.0) * du.Crossed(dv));
    const gp_Vec offset(rotation.Location(), point);
    const gp_Vec velocity     = -(vertical.Crossed(offset) + rise * vertical);
    const gp_Vec acceleration = vertical.Crossed(vertical.Crossed(offset));
    return BRepSolidSweep::ContactClearanceAcceleration(outward,
                                                        du,
                                                        dv,
                                                        duu,
                                                        duv,
                                                        dvv,
                                                        velocity,
                                                        acceleration);
  }

  gp_Pnt pose(const gp_Pnt& p, double t)
  {
    return p.Rotated(rotation, t).Translated((t * rise) * vertical);
  }

  Sample moving(const gp_Pnt& p, double t)
  {
    const gp_Pnt q = pose(p, t);
    return {q, vertical.Crossed(gp_Vec(rotation.Location(), q)) + rise * vertical};
  }

  template <class Function>
  std::vector<gp_Pnt> hermite(double end, int intervals, const Function& evaluate)
  {
    std::vector<gp_Pnt> poles;
    poles.reserve(3 * intervals + 1);
    const double step  = end / intervals;
    Sample       first = evaluate(0.0);
    poles.push_back(first.point);
    for (int i = 1; i <= intervals; ++i)
    {
      const Sample last = evaluate(i * step);
      poles.push_back(first.point.Translated((step / 3.0) * first.derivative));
      poles.push_back(last.point.Translated((-step / 3.0) * last.derivative));
      poles.push_back(last.point);
      first = last;
    }
    return poles;
  }

  void knots(NCollection_Array1<double>& values,
             NCollection_Array1<int>&    mults,
             double                      end,
             int                         intervals,
             int                         degree)
  {
    for (int i = 0; i <= intervals; ++i)
    {
      values(i + 1) = end * i / intervals;
      mults(i + 1)  = (i == 0 || i == intervals) ? degree + 1 : degree;
    }
  }

  template <class Column>
  occ::handle<Geom_BSplineSurface> surface(double        uEnd,
                                           int           uIntervals,
                                           int           uDegree,
                                           double        vEnd,
                                           int           vIntervals,
                                           int           vDegree,
                                           const Column& column)
  {
    const int                  uCount = uDegree * uIntervals + 1, vCount = vDegree * vIntervals + 1;
    NCollection_Array2<gp_Pnt> poles(1, uCount, 1, vCount);
    NCollection_Array2<double> weights(1, uCount, 1, vCount);
    for (int i = 0; i < uCount; ++i)
    {
      const auto values = column(i);
      for (int j = 0; j < vCount; ++j)
      {
        poles(i + 1, j + 1)   = values[j];
        weights(i + 1, j + 1) = uDegree == 2 && i == 1 ? std::sqrt(0.5) : 1.0;
      }
    }

    NCollection_Array1<double> uk(1, uIntervals + 1), vk(1, vIntervals + 1);
    NCollection_Array1<int>    um(1, uIntervals + 1), vm(1, vIntervals + 1);
    knots(uk, um, uEnd, uIntervals, uDegree);
    knots(vk, vm, vEnd, vIntervals, vDegree);
    return new Geom_BSplineSurface(poles, weights, uk, vk, um, vm, uDegree, vDegree);
  }

  occ::handle<Geom_SurfaceOfRevolution> endpointSurface(int part, double t)
  {
    const auto meridian = occ::down_cast<Geom_Curve>(bands[part].Curve->Copy());
    meridian->Rotate(rotation, t);
    meridian->Translate((t * rise) * vertical);
    return new Geom_SurfaceOfRevolution(meridian,
                                        gp_Ax1(pose(base, t), gp_Dir(axis.Rotated(rotation, t))));
  }

  gp_Pnt2d uv(int part, double u, bool last, bool upper)
  {
    const auto   e      = equation(profile(part, u));
    const double phase  = std::atan2(e.b, e.a);
    double       cosine = std::clamp(-e.c / std::hypot(e.a, e.b), -1.0, 1.0);
    if (std::abs(std::abs(cosine) - 1.0) < 1.0e-12)
      cosine = std::copysign(1.0, cosine);
    const double half = std::acos(cosine);
    const double azimuth =
      last ? phase + (upper ? half : -half) : phase + (upper ? 2.0 * pi - half : half);
    return gp_Pnt2d(azimuth, meridianParameter(part, u));
  }

  void addEndpointBand(BRepBuilderAPI_Sewing& sewing,
                       int                    part,
                       double                 a,
                       double                 b,
                       bool                   last,
                       bool                   hasContact)
  {
    const auto support = endpointSurface(part, last ? angle : 0.0);
    if (!hasContact)
    {
      const double            va = meridianParameter(part, a), vb = meridianParameter(part, b);
      BRepBuilderAPI_MakeFace face(support,
                                   0.0,
                                   2.0 * pi,
                                   std::min(va, vb),
                                   std::max(va, vb),
                                   1e-7);
      if (!face.IsDone())
        throw Standard_Failure("Could not build a full endpoint band");
      sewing.Add(face.Face());
      return;
    }

    constexpr int                  stations = 65;
    std::array<gp_Pnt2d, stations> lowerPoints, upperPoints;
    for (int i = 0; i < stations; ++i)
    {
      const double s     = double(i) / (stations - 1);
      const double u     = a + (b - a) * (0.5 * (1.0 - std::cos(pi * s)));
      auto         lower = uv(part, u, last, false);
      auto         upper = uv(part, u, last, true);
      if (i > 0)
      {
        // Lift both boundaries into the same continuous angular interval.
        const double shift =
          2.0 * pi * std::round((lowerPoints[i - 1].X() - lower.X()) / (2.0 * pi));
        lower.SetX(lower.X() + shift);
        upper.SetX(upper.X() + shift);
      }
      lowerPoints[i] = lower;
      upperPoints[i] = upper;
    }

    const auto makeBoundary = [&](const auto& samples, bool backwards) {
      occ::handle<NCollection_HArray1<gp_Pnt2d>> points =
        new NCollection_HArray1<gp_Pnt2d>(1, stations);
      for (int i = 0; i < stations; ++i)
        points->SetValue(i + 1, samples[backwards ? stations - 1 - i : i]);
      Geom2dAPI_Interpolate fit(points, false, 1e-10);
      fit.Perform();
      return fit.Curve();
    };

    BRepBuilderAPI_MakeWire outline;
    const auto              addStraight = [&](const gp_Pnt2d& first, const gp_Pnt2d& lastPoint) {
      if (first.Distance(lastPoint) > 1e-8)
      {
        GC_MakeSegment2d        segment(first, lastPoint);
        BRepBuilderAPI_MakeEdge edge(segment.Value(), support);
        outline.Add(edge.Edge());
      }
    };
    addStraight(lowerPoints.front(), upperPoints.front());
    BRepBuilderAPI_MakeEdge upper(makeBoundary(upperPoints, false), support);
    outline.Add(upper.Edge());
    addStraight(upperPoints.back(), lowerPoints.back());
    BRepBuilderAPI_MakeEdge lower(makeBoundary(lowerPoints, true), support);
    outline.Add(lower.Edge());
    if (!outline.IsDone())
      throw Standard_Failure("Could not join endpoint trim curves");
    TopoDS_Wire wire = outline.Wire();
    if (meridianParameter(part, b) < meridianParameter(part, a))
      wire.Reverse();
    BRepBuilderAPI_MakeFace face(support, wire, true);
    if (!face.IsDone() || !BRepLib::BuildCurves3d(face.Face(), 1e-5))
      throw Standard_Failure("Could not build a trimmed endpoint band");
    sewing.Add(face.Face());
  }

  void addDisk(BRepBuilderAPI_Sewing& sewing, bool last)
  {
    const double            t      = last ? angle : 0.0;
    const gp_Pnt            center = pose(base.Translated((last ? width : 0.0) * axis), t);
    const gp_Dir            normal(axis.Rotated(rotation, t));
    const gp_Dir            x(radialX.Rotated(rotation, t));
    BRepBuilderAPI_MakeEdge edge(
      gp_Circ(gp_Ax2(center, normal, x), last ? frontFaceRadius : backFaceRadius));
    BRepBuilderAPI_MakeWire wire(edge.Edge());
    BRepBuilderAPI_MakeFace face(wire.Wire());
    sewing.Add(face.Face());
  }

  double meridianRadius(int part) const
  {
    if (bands[part].Surface->GetType() == GeomAbs_Torus)
      return bands[part].Surface->Torus().MinorRadius();
    return std::numeric_limits<double>::infinity();
  }

  double meridianArcLength(int part, double a, double b) const
  {
    const double span = std::abs(meridianParameter(part, b) - meridianParameter(part, a));
    if (bands[part].Surface->GetType() == GeomAbs_Torus)
      return span * bands[part].Surface->Torus().MinorRadius();
    return span;
  }

  int contactMeridianIntervals(int part, double a, double b, int branch, double tolerance)
  {
    constexpr int aSamples = 16;
    gp_Pnt        aPoints[aSamples + 1];
    for (int i = 0; i <= aSamples; ++i)
      aPoints[i] = contact(part, a, b, branch, static_cast<double>(i) / aSamples);

    double aLength   = 0.0;
    double aTightest = std::numeric_limits<double>::infinity();
    for (int i = 1; i < aSamples; ++i)
    {
      const double aSide = aPoints[i - 1].Distance(aPoints[i]);
      const double bSide = aPoints[i].Distance(aPoints[i + 1]);
      const double cSide = aPoints[i - 1].Distance(aPoints[i + 1]);
      aLength += aSide;
      const gp_Vec aAcross(aPoints[i - 1], aPoints[i]);
      const gp_Vec bAcross(aPoints[i - 1], aPoints[i + 1]);
      const double aArea = aAcross.Crossed(bAcross).Magnitude();
      if (aArea <= 1.0e-18 || aSide <= Precision::Confusion() || bSide <= Precision::Confusion()
          || cSide <= Precision::Confusion())
        continue;
      aTightest = std::min(aTightest, aSide * bSide * cSide / (2.0 * aArea));
    }
    aLength += aPoints[aSamples - 1].Distance(aPoints[aSamples]);
    return BRepSolidSweep::ChordIntervals(aLength, aTightest, tolerance);
  }

  int screwIntervals(int part, double a, double b, int branch, double tolerance)
  {
    int count = 1;
    for (int i = 0; i <= 4; ++i)
    {
      const gp_Pnt point = contact(part, a, b, branch, i / 4.0);
      const gp_Vec from(rotation.Location(), point);
      const double axial  = from.Dot(vertical);
      const double rho    = (from - axial * vertical).Magnitude();
      const double speed  = std::hypot(rho, rise);
      const double arc    = std::abs(angle) * speed;
      const double radius = rho > Precision::Confusion() ? (rho * rho + rise * rise) / rho
                                                         : std::numeric_limits<double>::infinity();
      count               = std::max(count, BRepSolidSweep::ChordIntervals(arc, radius, tolerance));
    }
    return count;
  }

  void addInterval(BRepBuilderAPI_Sewing& sewing,
                   int                    part,
                   double                 a,
                   double                 b,
                   int&                   count,
                   double                 tolerance)
  {
    const bool useCurvature    = gScrewSchedule == BRepSolidSweep::ScrewStationSchedule::Curvature;
    const int  uniformMeridian = 32;
    const int  uniformScrew    = std::max(timeIntervals, 1);
    for (int branch : {-1, 1})
    {
      const auto pointAt         = [&](double u) { return contact(part, a, b, branch, u); };
      const int  surfaceMeridian = BRepSolidSweep::ChordIntervals(meridianArcLength(part, a, b),
                                                                  meridianRadius(part),
                                                                  tolerance);
      const int  contactMeridian = contactMeridianIntervals(part, a, b, branch, tolerance);
      // A toroidal corner spends points from its minor radius. A straight face
      // has infinite surface radius, so the contact curve supplies the radius.
      const int proposedMeridian = !useCurvature                         ? uniformMeridian
                                   : std::isfinite(meridianRadius(part)) ? surfaceMeridian
                                                                         : contactMeridian;
      const int proposedScrew =
        useCurvature ? screwIntervals(part, a, b, branch, tolerance) : uniformScrew;
      const int meridianCap = uniformMeridian * 8;
      const int screwCap    = uniformScrew * 8;
      int       stations    = std::clamp(proposedMeridian, 1, meridianCap);
      int       steps       = std::clamp(proposedScrew, 1, screwCap);
      bool      accepted    = false;
      for (int attempt = 0; attempt < 12 && !accepted; ++attempt)
      {
        const auto uPoles    = hermite(1.0, stations, [&](double u) {
          const double delta = 1e-4;
          const double lo = std::max(0.0, u - delta), hi = std::min(1.0, u + delta);
          return Sample{pointAt(u), gp_Vec(pointAt(lo), pointAt(hi)) / (hi - lo)};
        });
        const auto patch     = surface(1.0, stations, 3, angle, steps, 3, [&](int i) {
          return hermite(angle, steps, [&](double t) { return moving(uPoles[i], t); });
        });
        const auto deviation = [&](int uCount, int tCount, bool stopAtTolerance) {
          double worst = 0.0;
          for (int i = 0; i < uCount; ++i)
          {
            const double u     = (i + 0.5) / static_cast<double>(uCount);
            const gp_Pnt point = pointAt(u);
            for (int j = 0; j < tCount; ++j)
            {
              const double t = angle * (j + 0.5) / static_cast<double>(tCount);
              worst          = std::max(worst, pose(point, t).Distance(patch->Value(u, t)));
              if (stopAtTolerance && worst > tolerance)
                return worst;
            }
          }
          return worst;
        };
        // Mid-knot samples are the historical gate. The dense probe visits
        // points that grid does not, and is required for a curvature guess.
        const double coarse = deviation(3 * stations, steps, true);
        if (coarse > tolerance)
        {
          const int nextStations = std::min(stations * 2, meridianCap);
          const int nextSteps    = std::min(steps * 2, screwCap);
          if (nextStations == stations && nextSteps == steps)
            break;
          stations = nextStations;
          steps    = nextSteps;
          continue;
        }
        const double dense = deviation(8 * stations, 4 * steps, useCurvature);
        if (useCurvature && dense > tolerance)
        {
          const int nextStations = std::min(stations * 2, meridianCap);
          const int nextSteps    = std::min(steps * 2, screwCap);
          if (nextStations == stations && nextSteps == steps)
            break;
          stations = nextStations;
          steps    = nextSteps;
          continue;
        }

        BRepBuilderAPI_MakeFace face(patch, 1.0e-7);
        if (!face.IsDone() || !BRepCheck_Analyzer(face.Face()).IsValid())
          throw Standard_Failure("Could not build a valid screw envelope patch");
        sewing.Add(face.Face());
        ++count;
        gScrewFits.push_back(
          BRepSolidSweep::ScrewPatchFit{static_cast<int>(bands[part].Surface->GetType()),
                                        surfaceMeridian,
                                        contactMeridian,
                                        proposedMeridian,
                                        proposedScrew,
                                        uniformMeridian,
                                        uniformScrew,
                                        attempt + 1,
                                        stations,
                                        steps,
                                        coarse,
                                        dense});
        accepted = true;
      }
      if (!accepted)
        throw Standard_ConstructionError("Screw envelope missed the sampled fit tolerance");
    }
  }

  TopoDS_Solid build(double tolerance)
  {
    validateMeridian();
    if (std::abs(2.0 * pi * rise) <= width + 2.0 * maxRadius)
    {
      throw Standard_ConstructionError("Overlapping screw turns require envelope trimming");
    }
    const gp_Vec reference = vertical.Crossed(gp_Vec(rotation.Location(), base)) + rise * vertical;
    // The start exposes the back cap and the end exposes the front cap.
    // Other cap-sign cases require a different endpoint selection.
    if (axis.Dot(reference)
        <= maxRadius * axis.Crossed(vertical).Magnitude() + Precision::Confusion())
    {
      throw Standard_ConstructionError(
        "Screw sweep requires positive axial velocity over both caps");
    }

    BRepBuilderAPI_Sewing sewing(0.8 * tolerance);
    int                   count = 0;
    for (int part = 0; part < static_cast<int>(bands.size()); ++part)
    {
      std::vector<double> boundaries{0.0};
      for (int i = 1; i <= 128; ++i)
      {
        const double a = (i - 1) / 128.0, b = i / 128.0;
        if ((discriminant(part, a) >= 0.0) != (discriminant(part, b) >= 0.0))
          boundaries.push_back(root(part, a, b));
      }
      boundaries.push_back(1.0);
      for (std::size_t i = 1; i < boundaries.size(); ++i)
      {
        const double a = boundaries[i - 1], b = boundaries[i], mid = 0.5 * (a + b);
        if (discriminant(part, mid) > 0.0)
        {
          addInterval(sewing, part, a, b, count, tolerance);
          addEndpointBand(sewing, part, a, b, false, true);
          addEndpointBand(sewing, part, a, b, true, true);
        }
        else
        {
          addEndpointBand(sewing, part, a, b, equation(profile(part, mid)).c > 0.0, false);
        }
      }
    }
    if (count == 0)
    {
      throw Standard_ConstructionError("Screw envelope has no contact patches");
    }

    addDisk(sewing, false);
    addDisk(sewing, true);
    sewing.Perform();
    const TopoDS_Shape& shell = sewing.SewedShape();
    if (shell.IsNull() || shell.ShapeType() != TopAbs_SHELL || sewing.NbFreeEdges() != 0
        || sewing.NbMultipleEdges() != 0)
    {
      throw Standard_ConstructionError("Screw envelope did not close");
    }
    TopoDS_Solid result = BRepBuilderAPI_MakeSolid(TopoDS::Shell(shell)).Solid();
    if (!BRepLib::OrientClosedSolid(result) || !BRepCheck_Analyzer(result).IsValid())
    {
      throw Standard_ConstructionError("Screw envelope did not produce a valid solid");
    }
    return result;
  }
};

#include "BRepSolidSweep_ScrewTrim.pxx"
} // namespace

TopoDS_Solid BRepSolidSweep_BuildContactEnvelope(
  const TopoDS_Solid&                      theSolid,
  const BRepSolidSweep_Sweep::ScrewMotion& theMotion);

TopoDS_Solid BRepSolidSweep_Sweep::Sweep(const TopoDS_Solid& theSolid, const ScrewMotion& theMotion)
{
  BRepSolidSweep::ClearScrewPatchFits();
  if (!std::isfinite(theMotion.Angle) || theMotion.Angle < 0.0 || theMotion.Angle > 4.0 * pi
      || !std::isfinite(theMotion.AdvancePerTurn) || !std::isfinite(theMotion.SurfaceTolerance)
      || theMotion.SurfaceTolerance <= Precision::Confusion())
  {
    throw Standard_ConstructionError("Invalid screw motion or surface tolerance");
  }
  std::optional<AxisymmetricBody> aBody;
  try
  {
    aBody.emplace(readBody(theSolid, theMotion));
  }
  catch (const Standard_ConstructionError&)
  {
    return BRepSolidSweep_BuildContactEnvelope(theSolid, theMotion);
  }
  const auto& body = *aBody;
  if (theMotion.Angle <= Precision::Angular())
  {
    return theSolid;
  }
  ScrewEnvelope envelope(body, theMotion);
  const gp_Vec  aVelocity =
    envelope.vertical.Crossed(gp_Vec(envelope.rotation.Location(), envelope.base))
    + envelope.rise * envelope.vertical;
  const bool needsCapTrim =
    envelope.axis.Dot(aVelocity)
    <= envelope.maxRadius * envelope.axis.Crossed(envelope.vertical).Magnitude()
         + Precision::Confusion();
  if (needsCapTrim && std::abs(envelope.rise) > Precision::Confusion())
  {
    IntersectingScrewEnvelope trimmed(body, theMotion);
    return trimmed.buildTrimmed(theSolid, theMotion.SurfaceTolerance);
  }

  return envelope.build(theMotion.SurfaceTolerance);
}

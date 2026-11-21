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

#ifndef _BRepSolidSweep_MovingRibbons_HeaderFile
#define _BRepSolidSweep_MovingRibbons_HeaderFile

#include <GeomAPI_Interpolate.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <NCollection_HArray1.hxx>
#include <Standard_ConstructionError.hxx>

#include <cmath>
#include <map>
#include <vector>

namespace BRepSolidSweep
{
//! A contact curve at one station, with consistently ordered ends.
//! Equal junction indices mark ends that must follow the same trajectory.
struct MovingGenerator
{
  occ::handle<Geom_BoundedCurve> Curve;
  int                            FirstJunction;
  int                            LastJunction;
};

inline occ::handle<NCollection_HArray1<double>> RibbonParameters(int theIntervals)
{
  occ::handle<NCollection_HArray1<double>> aParams =
    new NCollection_HArray1<double>(1, theIntervals + 1);
  for (int i = 0; i <= theIntervals; ++i)
    aParams->SetValue(i + 1, double(i) / theIntervals);
  return aParams;
}

inline occ::handle<Geom_BSplineCurve> InterpolateRibbonCurve(
  const occ::handle<NCollection_HArray1<gp_Pnt>>& thePoints,
  const occ::handle<NCollection_HArray1<double>>& theParameters)
{
  GeomAPI_Interpolate anInterpolation(thePoints, theParameters, false, gp::Resolution());
  anInterpolation.Perform();
  if (!anInterpolation.IsDone())
    throw Standard_ConstructionError("Could not interpolate moving contact");
  return anInterpolation.Curve();
}

//! Fit continuous tensor-product patches from consistently identified contact
//! curves in world coordinates. Station is normalized to [0, 1]. Refine the
//! generator and station directions independently, checking between fit points.
template <class StationEvaluator>
std::vector<occ::handle<Geom_BSplineSurface>> FitMovingRibbons(StationEvaluator theEvaluate,
                                                               double           theTolerance)
{
  if (!(theTolerance > 0.0) || !std::isfinite(theTolerance))
    throw Standard_ConstructionError("Invalid moving contact tolerance");

  struct Station
  {
    std::vector<MovingGenerator> Generators;
    std::map<int, gp_Pnt>        Junctions;
  };

  std::map<double, Station> aStations;
  const auto                stationAt = [&](double theTime) -> const Station& {
    auto [anEntry, isNew] = aStations.try_emplace(theTime);
    auto& aStation        = anEntry->second;
    if (isNew)
    {
      aStation.Generators    = theEvaluate(theTime);
      const auto& aReference = aStations.at(0.0).Generators;
      if (aStation.Generators.empty() || aStation.Generators.size() != aReference.size())
        throw Standard_ConstructionError("Moving contact topology changed");

      const auto addJunction = [&](int theIndex, const gp_Pnt& thePoint) {
        auto [aJunction, isFirst] = aStation.Junctions.try_emplace(theIndex, thePoint);
        if (!isFirst && aJunction->second.Distance(thePoint) > 0.25 * theTolerance)
          throw Standard_ConstructionError("Moving contact junction separated");
      };
      for (size_t i = 0; i < aStation.Generators.size(); ++i)
      {
        const auto& aGenerator = aStation.Generators[i];
        if (aGenerator.Curve.IsNull() || aGenerator.FirstJunction != aReference[i].FirstJunction
            || aGenerator.LastJunction != aReference[i].LastJunction)
          throw Standard_ConstructionError("Moving contact correspondence changed");
        addJunction(aGenerator.FirstJunction, aGenerator.Curve->StartPoint());
        addJunction(aGenerator.LastJunction, aGenerator.Curve->EndPoint());
      }
    }
    return aStation;
  };

  const auto pointAt = [](const Station& theStation, size_t theIndex, double theU) {
    const auto& aGenerator = theStation.Generators[theIndex];
    if (theU == 0.0)
      return theStation.Junctions.at(aGenerator.FirstJunction);
    if (theU == 1.0)
      return theStation.Junctions.at(aGenerator.LastJunction);
    const auto& aCurve = aGenerator.Curve;
    return aCurve->Value(aCurve->FirstParameter()
                         + theU * (aCurve->LastParameter() - aCurve->FirstParameter()));
  };

  const auto fit = [&](size_t aGen, int aNu, int aNv) {
    const auto                                  aUParams = RibbonParameters(aNu);
    const auto                                  aVParams = RibbonParameters(aNv);
    std::vector<occ::handle<Geom_BSplineCurve>> aRows;
    occ::handle<NCollection_HArray1<gp_Pnt>>    aRow = new NCollection_HArray1<gp_Pnt>(1, aNu + 1);
    for (int j = 1; j <= aNv + 1; ++j)
    {
      const auto& aStation = stationAt(aVParams->Value(j));
      for (int i = 1; i <= aNu + 1; ++i)
        aRow->SetValue(i, pointAt(aStation, aGen, aUParams->Value(i)));
      aRows.push_back(InterpolateRibbonCurve(aRow, aUParams));
    }

    for (int j = 0; j <= aNv; ++j)
    {
      const auto& aStation = stationAt(aVParams->Value(j + 1));
      for (int i = 1; i < 4 * aNu; ++i)
      {
        const double u = double(i) / (4 * aNu);
        if (aRows[j]->Value(u).Distance(pointAt(aStation, aGen, u)) > 0.25 * theTolerance)
          return occ::handle<Geom_BSplineSurface>();
      }
    }

    std::vector<occ::handle<Geom_BSplineCurve>> aColumns;
    occ::handle<NCollection_HArray1<gp_Pnt>> aColumn = new NCollection_HArray1<gp_Pnt>(1, aNv + 1);
    for (int i = 1; i <= aRows.front()->NbPoles(); ++i)
    {
      for (int j = 1; j <= aNv + 1; ++j)
        aColumn->SetValue(j, aRows[j - 1]->Pole(i));
      aColumns.push_back(InterpolateRibbonCurve(aColumn, aVParams));
    }

    const auto&                aU = aRows.front();
    const auto&                aV = aColumns.front();
    NCollection_Array2<gp_Pnt> aPoles(1, aU->NbPoles(), 1, aV->NbPoles());
    for (int i = 1; i <= aU->NbPoles(); ++i)
      for (int j = 1; j <= aV->NbPoles(); ++j)
        aPoles.SetValue(i, j, aColumns[i - 1]->Pole(j));

    return occ::handle<Geom_BSplineSurface>(new Geom_BSplineSurface(aPoles,
                                                                    aU->Knots(),
                                                                    aV->Knots(),
                                                                    aU->Multiplicities(),
                                                                    aV->Multiplicities(),
                                                                    aU->Degree(),
                                                                    aV->Degree()));
  };

  const size_t  aCount        = stationAt(0.0).Generators.size();
  constexpr int aMaxIntervals = 128;
  // Sharing the station knots preserves identical junction trajectories.
  // Each patch can use its own generator knots.
  for (int aNv = 8; aNv <= aMaxIntervals; aNv *= 2)
  {
    std::vector<occ::handle<Geom_BSplineSurface>> aSurfaces;
    bool                                          needsMoreStations = false;
    for (size_t aGen = 0; aGen < aCount && !needsMoreStations; ++aGen)
    {
      occ::handle<Geom_BSplineSurface> aSurface;
      int                              aNu = 4;
      for (; aNu <= aMaxIntervals; aNu *= 2)
      {
        aSurface = fit(aGen, aNu, aNv);
        if (!aSurface.IsNull())
          break;
      }
      if (aSurface.IsNull())
        throw Standard_ConstructionError("Moving generator interpolation exceeds tolerance");

      for (int j = 1; j < 4 * aNv && !needsMoreStations; ++j)
      {
        if (j % 4 == 0)
          continue;
        const double v        = double(j) / (4 * aNv);
        const auto&  aStation = stationAt(v);
        for (int i = 0; i <= 4 * aNu; ++i)
        {
          const double u = double(i) / (4 * aNu);
          if (aSurface->Value(u, v).Distance(pointAt(aStation, aGen, u)) > 0.5 * theTolerance)
          {
            needsMoreStations = true;
            break;
          }
        }
      }
      aSurfaces.push_back(aSurface);
    }
    if (!needsMoreStations)
      return aSurfaces;
  }
  throw Standard_ConstructionError("Moving contact interpolation exceeds tolerance");
}
} // namespace BRepSolidSweep

#endif

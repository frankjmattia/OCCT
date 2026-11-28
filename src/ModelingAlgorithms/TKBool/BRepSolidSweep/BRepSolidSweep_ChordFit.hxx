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

#ifndef _BRepSolidSweep_ChordFit_HeaderFile
#define _BRepSolidSweep_ChordFit_HeaderFile

#include <Standard_Macro.hxx>

#include <cmath>
#include <limits>
#include <vector>

namespace BRepSolidSweep
{
//! Intervals of a circular arc whose chord stays within theTolerance.
//! A non-positive or infinite radius is a straight span: the two endpoints.
inline int ChordIntervals(double theArcLength, double theRadius, double theTolerance)
{
  if (!(theArcLength > 0.0) || !std::isfinite(theArcLength))
    return 1;

  if (!(theRadius > 0.0) || !std::isfinite(theRadius) || !(theTolerance > 0.0))
    return 1;

  const double aDisc = 2.0 * theRadius * theTolerance - theTolerance * theTolerance;
  if (!(aDisc > 0.0))
    return 1;

  const double aChord = 2.0 * std::sqrt(aDisc);
  if (!(aChord > 0.0) || !std::isfinite(aChord))
    return 1;
  const double aCount = std::ceil(theArcLength / aChord);
  if (aCount < 1.0)
    return 1;
  if (aCount > 1024.0)
    return 1024;
  return static_cast<int>(aCount);
}

//! Uniform repeats the historical fixed station schedule.
//! Curvature sets the first guess from ChordIntervals and still doubles on failure.
enum class ScrewStationSchedule
{
  Uniform,
  Curvature
};

//! One accepted contact patch from the screw envelope fitter.
struct ScrewPatchFit
{
  int    SurfaceType      = 0;
  int    SurfaceMeridian  = 0;
  int    ContactMeridian  = 0;
  int    ProposedMeridian = 0;
  int    ProposedScrew    = 0;
  int    UniformMeridian  = 0;
  int    UniformScrew     = 0;
  int    Attempts         = 0;
  int    AcceptedMeridian = 0;
  int    AcceptedScrew    = 0;
  double CoarseDeviation  = 0.0;
  double DenseDeviation   = 0.0;
};

Standard_EXPORT void                 SetScrewStationSchedule(ScrewStationSchedule theSchedule);
Standard_EXPORT ScrewStationSchedule ActiveScrewStationSchedule();
Standard_EXPORT void                 ClearScrewPatchFits();
Standard_EXPORT const std::vector<ScrewPatchFit>& ScrewPatchFits();
} // namespace BRepSolidSweep

#endif

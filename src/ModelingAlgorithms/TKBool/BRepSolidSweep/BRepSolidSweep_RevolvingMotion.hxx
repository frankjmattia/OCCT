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

#ifndef _BRepSolidSweep_RevolvingMotion_HeaderFile
#define _BRepSolidSweep_RevolvingMotion_HeaderFile

#include <Precision.hxx>
#include <Standard_ConstructionError.hxx>
#include <gp_Ax1.hxx>
#include <gp_Trsf.hxx>
#include <cmath>
#include <optional>

namespace BRepSolidSweep
{
//! Rigid motion consisting of a linear drift in the rotating frame. The
//! angular rate can be constant or inversely proportional to a linear law.
//! Its parameter is measured from the initial pose (identity).
class RevolvingMotion
{
public:
  RevolvingMotion(const gp_Ax1& theAxis,
                  const gp_Vec& theDrift,
                  double        theInitialRate,
                  double        theRateDenominatorSlope)
      : myAxis(theAxis),
        myDrift(theDrift),
        myInitialRate(theInitialRate),
        myRateDenominatorSlope(theRateDenominatorSlope)
  {
    if (!std::isfinite(myInitialRate) || !std::isfinite(myRateDenominatorSlope)
        || !std::isfinite(myDrift.SquareMagnitude()))
      throw Standard_ConstructionError("Non-finite revolving motion");
  }

  double Phase(double theStation) const
  {
    const double aRelativeChange = checkedRelativeChange(theStation);
    if (aRelativeChange == 0.0)
      return myInitialRate * theStation;
    return myInitialRate * theStation * (std::log1p(aRelativeChange) / aRelativeChange);
  }

  double AngularRate(double theStation) const
  {
    return myInitialRate / (1.0 + checkedRelativeChange(theStation));
  }

  gp_Trsf Placement(double theStation) const
  {
    gp_Trsf aRotation, aTranslation;
    aRotation.SetRotation(myAxis, Phase(theStation));
    aTranslation.SetTranslation(theStation * myDrift);
    return aRotation * aTranslation;
  }

  //! The angular and linear terms of the instantaneous velocity field.
  //! Velocity at a world point p is Angular x (p - AxisLocation) + Linear.
  struct VelocityField
  {
    gp_Pnt AxisLocation;
    gp_Vec Angular;
    gp_Vec Linear;

    gp_Vec At(const gp_Pnt& thePoint) const
    {
      return Angular.Crossed(gp_Vec(AxisLocation, thePoint)) + Linear;
    }
  };

  VelocityField Velocity(double theStation) const
  {
    return {myAxis.Location(),
            AngularRate(theStation) * gp_Vec(myAxis.Direction()),
            myDrift.Rotated(myAxis, Phase(theStation))};
  }

  struct ReferenceScrew
  {
    gp_Ax1 Axis;
    double Rise;
  };

  //! Velocity of the moving solid pulled back into its initial frame.
  VelocityField ReferenceVelocity(double theStation) const
  {
    return {myAxis.Location().Translated(-theStation * myDrift),
            AngularRate(theStation) * gp_Vec(myAxis.Direction()),
            myDrift};
  }

  //! Instantaneous screw in the initial frame, with a positive angular
  //! speed about its returned axis. No axis is returned for translation.
  std::optional<ReferenceScrew> ScrewAt(double theStation) const
  {
    const auto   aField = ReferenceVelocity(theStation);
    const double aSpeed = aField.Angular.Magnitude();
    if (aSpeed <= Precision::Angular())
      return std::nullopt;

    const gp_Dir aDirection(aField.Angular);
    const gp_Vec aUnit(aDirection);
    const gp_Pnt aCenter = aField.AxisLocation.Translated(aUnit.Crossed(aField.Linear) / aSpeed);
    return ReferenceScrew{gp_Ax1(aCenter, aDirection), aField.Linear.Dot(aUnit) / aSpeed};
  }

private:
  double checkedRelativeChange(double theStation) const
  {
    const double aRelativeChange = myRateDenominatorSlope * theStation;
    if (!std::isfinite(aRelativeChange) || aRelativeChange <= -1.0)
      throw Standard_ConstructionError("Station outside revolving motion domain");
    return aRelativeChange;
  }

  gp_Ax1 myAxis;
  gp_Vec myDrift;
  double myInitialRate;
  double myRateDenominatorSlope;
};
} // namespace BRepSolidSweep

#endif

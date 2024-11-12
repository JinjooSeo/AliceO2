// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsTRD/Constants.h"
#include "TRDBase/TrackletTransformer.h"
#include "TMath.h"

using namespace o2::trd;
using namespace o2::trd::constants;

void TrackletTransformer::init()
{
  mGeo = Geometry::instance();
  mGeo->createPadPlaneArray();
  mGeo->createClusterMatrixArray();

  // 3.35 cm
  mXAnode = mGeo->cdrHght() + mGeo->camHght() / 2;
}

float TrackletTransformer::calculateZ(int padrow, const PadPlane* padPlane) const
{
  double rowPos = padPlane->getRowPos(padrow);
  double rowSize = padPlane->getRowSize(padrow);
  double middleRowPos = padPlane->getRowPos(padPlane->getNrows() / 2);

  return rowPos - rowSize / 2. - middleRowPos;
}

float TrackletTransformer::calculateDy(float vDrift, float exb, int slope, const PadPlane* padPlane) const
{
  double padWidth = padPlane->getWidthIPad();

  // dy = slope * nTimeBins * padWidth * GRANULARITYTRKLSLOPE;
  // nTimeBins should be number of timebins in drift region. 1 timebin is 100 nanosecond
  double rawDy = slope * ((mGeo->cdrHght() / vDrift) * 10.) * padWidth * GRANULARITYTRKLSLOPE / ADDBITSHIFTSLOPE;

  // NOTE: confirm the direction in which vDrift is measured/determined. Is it in x or in direction of drift?
  double lorentzCorrection = TMath::Tan(exb) * mXAnode;

  // assuming angle in Bailhache, fig. 4.17 would be positive in our calibration code
  double calibratedDy = rawDy - lorentzCorrection;

  return calibratedDy;
}

float TrackletTransformer::calculateX(int detector, float t0, float vDrift) const
{
  // The drift time can be estimated according to the expression:
  // t_{drift} = t_{bin} - t_{0} - t_{cause}(x) - t_{TC}(q_{i-1}, q_{i-2}, ...)
  // t_0 is the delay of the trigger signal
  // t_{cause}(x) is the causality delay between ionisation electrons hitting the anode and the registration of maximum signal by the electronics - it is due to the rising time of the TRF 
  // t_TC is the residual charge from previous (in time) clusters due to residual tails after tail cancellation. 
  
  // 5mm below cathode plane to reduce error propogation from tracklet fit and driftV
  // Inverse the x position to get the correct sign
  // 0.5 - mGeo->cdrHght() = 0.5 - 3.35 = -2.85
  // 10 is a sampling frequency
  float td = getTimebin(detector, 0.5 - mGeo->cdrHght())/10.;

  // correction for t0
  td -= t0;
  if(td < .2) return 0.;
  
  // TRF rising time (fitted) 
  // It should be absorbed by the t0. For the moment t0 is 0 for simulations.
  // t_{cause}(x) = 0.189
  // A.Bercuci (Mar 26 2009)
  td -= 0.189;

  // The drift velocity varies with the drift length (and distance to anode wire) as described by cell structure simulation. 
  // Thus one, in principle, can calculate iteratively the drift length from the expression:
  // x = t_{drift}(x)*v_{drift}(x)
  return td * vDrift;;
}

std::array<float, 3> TrackletTransformer::transformL2T(int detector, std::array<double, 3> point) const
{
  auto transformationMatrix = mGeo->getMatrixT2L(detector);

  ROOT::Math::Impl::Transform3D<double>::Point localPoint(point[0], point[1], point[2]);
  auto gobalPoint = transformationMatrix ^ localPoint;

  return {(float)gobalPoint.x(), (float)gobalPoint.y(), (float)gobalPoint.z()};
}

CalibratedTracklet TrackletTransformer::transformTracklet(Tracklet64 tracklet, bool trackingFrame) const
{
  auto detector = tracklet.getDetector();
  auto hcid = tracklet.getHCID();
  auto padrow = tracklet.getPadRow();
  auto column = tracklet.getColumn();
  int position;
  int slope;
  
  if (mApplyXOR) {
    position = tracklet.getPosition() ^ 0x80;
    if (position & (1 << (constants::NBITSTRKLPOS - 1))) {
      position = -((~(position - 1)) & ((1 << constants::NBITSTRKLPOS) - 1));
    }
    slope = tracklet.getSlope() ^ 0x80;
    if (slope & (1 << (constants::NBITSTRKLSLOPE - 1))) {
      slope = -((~(slope - 1)) & ((1 << constants::NBITSTRKLSLOPE) - 1));
    }
  } else {
    position = tracklet.getPositionBinSigned();
    slope = tracklet.getSlopeBinSigned();
  }

  // calculate raw local chamber space point
  const auto padPlane = mGeo->getPadPlane(detector);

  // get calibration values
  float t0 = mCalT0->getT0(detector);
  float vDrift = mCalVdriftExB->getVdrift(detector);
  float exb = mCalVdriftExB->getExB(detector);

  // get position
  float x = calculateX(detector, t0, vDrift);
  float y = tracklet.getUncalibratedY(mApplyShift);
  float z = calculateZ(padrow, padPlane);

  // get slope
  float dy = calculateDy(vDrift, exb, slope, padPlane);
  
  // // apply ExB correction (https://github.com/alisw/AliRoot/blob/3bf76e5e3a8d88e1f7291637ee27a1c64327a9dc/TRD/TRDbase/AliTRDtransform.cxx#L307C5-L307C59)
  float calibratedY = y - (x * exb); 
  
  if (trackingFrame) {
    std::array<float, 3> sectorSpacePoint = transformL2T(detector, std::array<double, 3>{x, calibratedY, z});
    LOG(debug) << "x: " << sectorSpacePoint[0] << " | "
               << "y: " << sectorSpacePoint[1] << " | "
               << "z: " << sectorSpacePoint[2];
    return CalibratedTracklet(sectorSpacePoint[0], sectorSpacePoint[1], sectorSpacePoint[2], dy);
  } else {
    return CalibratedTracklet(x, calibratedY, z, dy); // local frame
  }
}

double TrackletTransformer::getTimebin(int detector, double x) const
{
  // calculate timebin from x position within chamber
  float vDrift = mCalVdriftExB->getVdrift(detector);
  double t0 = 4.0; // time (in timebins) of start of drift region

  double timebin;
  // x = 0 at anode plane and points toward pad plane.
  if (x < -mGeo->camHght() / 2) {
    // drift region
    timebin = t0 - (x + mGeo->camHght() / 2) / (vDrift * 0.1);
  } else {
    // anode region: very rough guess
    timebin = t0 - 1.0 + fabs(x);
  }

  return timebin;
}

namespace std {}
using namespace std;
#include "recovered_headersProjectHeaders.h"

#include "recovered_headersLinkDef.h"

#include "recovered_headersProjectDict.cxx"

struct DeleteObjectFunctor {
   template <typename T>
   void operator()(const T *ptr) const {
      delete ptr;
   }
   template <typename T, typename Q>
   void operator()(const std::pair<T,Q> &) const {
      // Do nothing
   }
   template <typename T, typename Q>
   void operator()(const std::pair<T,Q*> &ptr) const {
      delete ptr.second;
   }
   template <typename T, typename Q>
   void operator()(const std::pair<T*,Q> &ptr) const {
      delete ptr.first;
   }
   template <typename T, typename Q>
   void operator()(const std::pair<T*,Q*> &ptr) const {
      delete ptr.first;
      delete ptr.second;
   }
};

#ifndef tpcData_cxx
#define tpcData_cxx
tpcData::tpcData() {
}
tpcData &tpcData::operator=(const tpcData & rhs)
{
   // This is NOT a copy operator=. This is actually a move operator= (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   val = (const_cast<tpcData &>( rhs ).val);
   channel = (const_cast<tpcData &>( rhs ).channel);
   row = (const_cast<tpcData &>( rhs ).row);
   column = (const_cast<tpcData &>( rhs ).column);
   timestamp = (const_cast<tpcData &>( rhs ).timestamp);
   pedestal = (const_cast<tpcData &>( rhs ).pedestal);
   peddev = (const_cast<tpcData &>( rhs ).peddev);
   return *this;
}
tpcData::tpcData(const tpcData & rhs)
   : val(const_cast<tpcData &>( rhs ).val)
   , channel(const_cast<tpcData &>( rhs ).channel)
   , row(const_cast<tpcData &>( rhs ).row)
   , column(const_cast<tpcData &>( rhs ).column)
   , timestamp(const_cast<tpcData &>( rhs ).timestamp)
   , pedestal(const_cast<tpcData &>( rhs ).pedestal)
   , peddev(const_cast<tpcData &>( rhs ).peddev)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
}
tpcData::~tpcData() {
}
#endif // tpcData_cxx

#ifndef CentroidData_cxx
#define CentroidData_cxx
CentroidData::CentroidData() {
}
CentroidData &CentroidData::operator=(const CentroidData & rhs)
{
   // This is NOT a copy operator=. This is actually a move operator= (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   y = (const_cast<CentroidData &>( rhs ).y);
   x = (const_cast<CentroidData &>( rhs ).x);
   z = (const_cast<CentroidData &>( rhs ).z);
   integrated_ADC_amplitude = (const_cast<CentroidData &>( rhs ).integrated_ADC_amplitude);
   sigma_x = (const_cast<CentroidData &>( rhs ).sigma_x);
   sigma_y = (const_cast<CentroidData &>( rhs ).sigma_y);
   sigma_z = (const_cast<CentroidData &>( rhs ).sigma_z);
   used_in_fit = (const_cast<CentroidData &>( rhs ).used_in_fit);
   member_tpc_idx = (const_cast<CentroidData &>( rhs ).member_tpc_idx);
   CentroidData &modrhs = const_cast<CentroidData &>( rhs );
   modrhs.member_tpc_idx.clear();
   return *this;
}
CentroidData::CentroidData(const CentroidData & rhs)
   : y(const_cast<CentroidData &>( rhs ).y)
   , x(const_cast<CentroidData &>( rhs ).x)
   , z(const_cast<CentroidData &>( rhs ).z)
   , integrated_ADC_amplitude(const_cast<CentroidData &>( rhs ).integrated_ADC_amplitude)
   , sigma_x(const_cast<CentroidData &>( rhs ).sigma_x)
   , sigma_y(const_cast<CentroidData &>( rhs ).sigma_y)
   , sigma_z(const_cast<CentroidData &>( rhs ).sigma_z)
   , used_in_fit(const_cast<CentroidData &>( rhs ).used_in_fit)
   , member_tpc_idx(const_cast<CentroidData &>( rhs ).member_tpc_idx)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   CentroidData &modrhs = const_cast<CentroidData &>( rhs );
   modrhs.member_tpc_idx.clear();
}
CentroidData::~CentroidData() {
}
#endif // CentroidData_cxx

#ifndef TrackData_cxx
#define TrackData_cxx
TrackData::TrackData() {
}
TrackData &TrackData::operator=(const TrackData & rhs)
{
   // This is NOT a copy operator=. This is actually a move operator= (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   nPoints = (const_cast<TrackData &>( rhs ).nPoints);
   slope_xy = (const_cast<TrackData &>( rhs ).slope_xy);
   intercept_xy = (const_cast<TrackData &>( rhs ).intercept_xy);
   slope_zy = (const_cast<TrackData &>( rhs ).slope_zy);
   intercept_zy = (const_cast<TrackData &>( rhs ).intercept_zy);
   y = (const_cast<TrackData &>( rhs ).y);
   x = (const_cast<TrackData &>( rhs ).x);
   z = (const_cast<TrackData &>( rhs ).z);
   charge = (const_cast<TrackData &>( rhs ).charge);
   sigmas_x = (const_cast<TrackData &>( rhs ).sigmas_x);
   sigmas_z = (const_cast<TrackData &>( rhs ).sigmas_z);
   sigmas_y = (const_cast<TrackData &>( rhs ).sigmas_y);
   slope_xy_err = (const_cast<TrackData &>( rhs ).slope_xy_err);
   intercept_xy_err = (const_cast<TrackData &>( rhs ).intercept_xy_err);
   slope_zy_err = (const_cast<TrackData &>( rhs ).slope_zy_err);
   intercept_zy_err = (const_cast<TrackData &>( rhs ).intercept_zy_err);
   cov_xy_is = (const_cast<TrackData &>( rhs ).cov_xy_is);
   cov_zy_is = (const_cast<TrackData &>( rhs ).cov_zy_is);
   chi2_xy = (const_cast<TrackData &>( rhs ).chi2_xy);
   ndf_xy = (const_cast<TrackData &>( rhs ).ndf_xy);
   chi2ndf_xy = (const_cast<TrackData &>( rhs ).chi2ndf_xy);
   p_xy = (const_cast<TrackData &>( rhs ).p_xy);
   chi2_zy = (const_cast<TrackData &>( rhs ).chi2_zy);
   ndf_zy = (const_cast<TrackData &>( rhs ).ndf_zy);
   chi2ndf_zy = (const_cast<TrackData &>( rhs ).chi2ndf_zy);
   p_zy = (const_cast<TrackData &>( rhs ).p_zy);
   TrackData &modrhs = const_cast<TrackData &>( rhs );
   modrhs.y.clear();
   modrhs.x.clear();
   modrhs.z.clear();
   modrhs.charge.clear();
   modrhs.sigmas_x.clear();
   modrhs.sigmas_z.clear();
   modrhs.sigmas_y.clear();
   return *this;
}
TrackData::TrackData(const TrackData & rhs)
   : nPoints(const_cast<TrackData &>( rhs ).nPoints)
   , slope_xy(const_cast<TrackData &>( rhs ).slope_xy)
   , intercept_xy(const_cast<TrackData &>( rhs ).intercept_xy)
   , slope_zy(const_cast<TrackData &>( rhs ).slope_zy)
   , intercept_zy(const_cast<TrackData &>( rhs ).intercept_zy)
   , y(const_cast<TrackData &>( rhs ).y)
   , x(const_cast<TrackData &>( rhs ).x)
   , z(const_cast<TrackData &>( rhs ).z)
   , charge(const_cast<TrackData &>( rhs ).charge)
   , sigmas_x(const_cast<TrackData &>( rhs ).sigmas_x)
   , sigmas_z(const_cast<TrackData &>( rhs ).sigmas_z)
   , sigmas_y(const_cast<TrackData &>( rhs ).sigmas_y)
   , slope_xy_err(const_cast<TrackData &>( rhs ).slope_xy_err)
   , intercept_xy_err(const_cast<TrackData &>( rhs ).intercept_xy_err)
   , slope_zy_err(const_cast<TrackData &>( rhs ).slope_zy_err)
   , intercept_zy_err(const_cast<TrackData &>( rhs ).intercept_zy_err)
   , cov_xy_is(const_cast<TrackData &>( rhs ).cov_xy_is)
   , cov_zy_is(const_cast<TrackData &>( rhs ).cov_zy_is)
   , chi2_xy(const_cast<TrackData &>( rhs ).chi2_xy)
   , ndf_xy(const_cast<TrackData &>( rhs ).ndf_xy)
   , chi2ndf_xy(const_cast<TrackData &>( rhs ).chi2ndf_xy)
   , p_xy(const_cast<TrackData &>( rhs ).p_xy)
   , chi2_zy(const_cast<TrackData &>( rhs ).chi2_zy)
   , ndf_zy(const_cast<TrackData &>( rhs ).ndf_zy)
   , chi2ndf_zy(const_cast<TrackData &>( rhs ).chi2ndf_zy)
   , p_zy(const_cast<TrackData &>( rhs ).p_zy)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   TrackData &modrhs = const_cast<TrackData &>( rhs );
   modrhs.y.clear();
   modrhs.x.clear();
   modrhs.z.clear();
   modrhs.charge.clear();
   modrhs.sigmas_x.clear();
   modrhs.sigmas_z.clear();
   modrhs.sigmas_y.clear();
}
TrackData::~TrackData() {
}
#endif // TrackData_cxx


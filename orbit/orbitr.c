/* N3EMO Orbit Simulator routines  v3.7 */

/* Copyright (c) 1986,1987,1988,1989,1990 Robert W. Berger N3EMO
   May be freely distributed, provided this notice remains intact. */

#include <math.h>
#include <stdio.h>
#include <string.h>

#define SSPELLIPSE 0
/* If non-zero, use ellipsoidal earth model when calculating longitude,
 * latitude, and height.
 */

#ifndef PI
#define PI 3.1415926535897932384626433832795028841971693993751058
#endif

#ifdef PI2
#undef PI2
#endif

typedef double mat3x3[3][3];

#define PI2 (PI * 2)
#define MinutesPerDay (24 * 60.0)
#define SecondsPerDay (60 * MinutesPerDay)
#define HalfSecond (0.5 / SecondsPerDay)
#define EarthRadius 6378.16 /* Kilometers, at equator */

#define EarthFlat (1 / 298.25)         /* Earth Flattening Coeff. */
#define SiderealSolar 1.0027379093
#define SidRate (PI2 * SiderealSolar / SecondsPerDay) /* radians/second */
#define GM 398600                                 /* Kilometers^3/seconds^2 */
#define DegreesPerRadian (180 / PI)
#define RadiansPerDegree (PI / 180)
#define ABS(x) ((x) < 0 ? (-(x)) : (x))
#define SQR(x) ((x) * (x))

#define Epsilon (RadiansPerDegree / 3600) /* 1 arc second */
#define SunRadius 695000
#define SunSemiMajorAxis 149598845.0 /* Kilometers */

/* Sidereal time state. */
static double SidDay;
static double SidReference;

/* Keplerian elements for the sun. */
static double SunEpochTime;
static double SunInclination;
static double SunRAAN;
static double SunEccentricity;
static double SunArgPerigee;
static double SunMeanAnomaly;
static double SunMeanMotion;

/* Values for shadow geometry. */
static double SinPenumbra;
static double CosPenumbra;

/* Local helpers. */
static void GetTopocentric(double sat_x, double sat_y, double sat_z,
                           double site_x, double site_y, double site_z,
                           mat3x3 site_matrix, double *x, double *y,
                           double *z);
static void SPrintDate(char *str, long day_num);
static void SPrintDayOfWeek(char *str, long day_num);
static void SPrintTime(char *str, double time);
static void GetDate(long day_num, int *yearp, int *monthp, int *dayp);

long calls = 0;
long iters = 0;

/**
 * dumpstats - print Kepler() iteration statistics
 *
 * This is a debugging helper.
 */
void dumpstats(void)
{
   if (calls == 0) {
      printf("Average iterations = (no calls)\n");
      return;
   }

   printf("Average iterations = %lf\n", ((double)iters) / calls);
}

/**
 * Kepler - solve Kepler's equation
 * @mean_anomaly: time since last perigee in radians (PI2 = one orbit)
 * @eccentricity: eccentricity of orbit's ellipse
 *
 * Return: true anomaly in the range [0, PI2).
 */
double Kepler(double mean_anomaly, double eccentricity)
{
   double ecc_anom; /* Eccentric anomaly */
   double error;
   double true_anomaly;

   calls++;

   /* Initial guess. */
   ecc_anom = mean_anomaly;

   do {
      error = (ecc_anom - eccentricity * sin(ecc_anom) - mean_anomaly) /
              (1 - eccentricity * cos(ecc_anom));
      ecc_anom -= error;
      iters++;
   } while (ABS(error) >= Epsilon);

   if (ABS(ecc_anom - PI) < Epsilon)
      true_anomaly = PI;
   else
      true_anomaly = 2 * atan(sqrt((1 + eccentricity) / (1 - eccentricity)) *
                              tan(ecc_anom / 2));

   if (true_anomaly < 0)
      true_anomaly += PI2;

   return true_anomaly;
}

/**
 * GetSubSatPoint - compute sub-satellite point
 * @sat_x: satellite x coordinate (km)
 * @sat_y: satellite y coordinate (km)
 * @sat_z: satellite z coordinate (km)
 * @time: current time in days
 * @latitude: output latitude (radians)
 * @longitude: output longitude (radians), wrapped to [0, PI2)
 * @height: output height above EarthRadius (km)
 */
void GetSubSatPoint(double sat_x, double sat_y, double sat_z, double time,
                    double *latitude, double *longitude, double *height)
{
   double r;
   long i;

   r = sqrt(SQR(sat_x) + SQR(sat_y) + SQR(sat_z));

   *longitude = PI2 * ((time - SidDay) * SiderealSolar + SidReference) -
                atan2(sat_y, sat_x);

   /* i = floor(longitude / (2*pi)) */
   i = (long)(*longitude / PI2);
   if (i < 0)
      i--;

   *longitude -= (double)i * PI2;
   *latitude = atan(sat_z / sqrt(SQR(sat_x) + SQR(sat_y)));

#if SSPELLIPSE
   (void)height;
#else
   *height = r - EarthRadius;
#endif
}

/**
 * GetPrecession - compute orbital precession rates
 * @semi_major_axis: semi-major axis (km)
 * @eccentricity: orbital eccentricity
 * @inclination: orbital inclination (radians)
 * @raan_precession: output precession of RAAN (radians/day)
 * @perigee_precession: output precession of perigee (radians/day)
 */
void GetPrecession(double semi_major_axis, double eccentricity,
                   double inclination, double *raan_precession,
                   double *perigee_precession)
{
   *raan_precession =
      9.95 * pow(EarthRadius / semi_major_axis, 3.5) * cos(inclination) /
      SQR(1 - SQR(eccentricity)) * RadiansPerDegree;

   *perigee_precession =
      4.97 * pow(EarthRadius / semi_major_axis, 3.5) *
      (5 * SQR(cos(inclination)) - 1) / SQR(1 - SQR(eccentricity)) *
      RadiansPerDegree;
}

/**
 * GetSatPosition - compute satellite position and velocity
 * @epoch_time: epoch time (days)
 * @epoch_raan: RAAN at epoch (radians)
 * @epoch_arg_perigee: argument of perigee at epoch (radians)
 * @semi_major_axis: semi-major axis (km)
 * @inclination: orbital inclination (radians)
 * @eccentricity: orbital eccentricity
 * @raan_precession: RAAN precession rate
 * @perigee_precession: perigee precession rate
 * @time: time at which to compute the state (days)
 * @true_anomaly: true anomaly at @time (radians)
 * @x: output x coordinate (km)
 * @y: output y coordinate (km)
 * @z: output z coordinate (km)
 * @radius: output radius from geocenter (km)
 * @vx: output x velocity (km/s)
 * @vy: output y velocity (km/s)
 * @vz: output z velocity (km/s)
 */
void GetSatPosition(double epoch_time, double epoch_raan,
                    double epoch_arg_perigee, double semi_major_axis,
                    double inclination, double eccentricity,
                    double raan_precession, double perigee_precession,
                    double time, double true_anomaly, double *x, double *y,
                    double *z, double *radius, double *vx, double *vy,
                    double *vz)
{
   double raan;
   double arg_perigee;
   double xw, yw, vxw, vyw; /* In orbital plane */
   double tmp;
   double px, qx, py, qy, pz, qz; /* Escobal transformation 31 */
   double cos_arg_perigee, sin_arg_perigee;
   double cos_raan, sin_raan, cos_inclination, sin_inclination;

   *radius = semi_major_axis * (1 - SQR(eccentricity)) /
             (1 + eccentricity * cos(true_anomaly));

   xw = *radius * cos(true_anomaly);
   yw = *radius * sin(true_anomaly);

   tmp = sqrt(GM / (semi_major_axis * (1 - SQR(eccentricity))));
   vxw = -tmp * sin(true_anomaly);
   vyw = tmp * (cos(true_anomaly) + eccentricity);

   arg_perigee = epoch_arg_perigee + (time - epoch_time) * perigee_precession;
   raan = epoch_raan - (time - epoch_time) * raan_precession;

   cos_raan = cos(raan);
   sin_raan = sin(raan);
   cos_arg_perigee = cos(arg_perigee);
   sin_arg_perigee = sin(arg_perigee);
   cos_inclination = cos(inclination);
   sin_inclination = sin(inclination);

   px = cos_arg_perigee * cos_raan -
        sin_arg_perigee * sin_raan * cos_inclination;
   py = cos_arg_perigee * sin_raan +
        sin_arg_perigee * cos_raan * cos_inclination;
   pz = sin_arg_perigee * sin_inclination;
   qx = -sin_arg_perigee * cos_raan -
        cos_arg_perigee * sin_raan * cos_inclination;
   qy = -sin_arg_perigee * sin_raan +
        cos_arg_perigee * cos_raan * cos_inclination;
   qz = cos_arg_perigee * sin_inclination;

   /* Escobal, transformation #31. */
   *x = px * xw + qx * yw;
   *y = py * xw + qy * yw;
   *z = pz * xw + qz * yw;

   *vx = px * vxw + qx * vyw;
   *vy = py * vxw + qy * vyw;
   *vz = pz * vxw + qz * vyw;
}

/**
 * GetSitPosition - compute site position/velocity and transform matrix
 * @site_lat: site latitude (radians)
 * @site_long: site longitude (radians)
 * @site_elevation: site elevation above sea level (km)
 * @current_time: current time (days)
 * @site_x: output site x coordinate (km)
 * @site_y: output site y coordinate (km)
 * @site_z: output site z coordinate (km)
 * @site_vx: output site x velocity (km/s)
 * @site_vy: output site y velocity (km/s)
 * @site_matrix: output matrix converting geocentric to topocentric coords
 */
void GetSitPosition(double site_lat, double site_long, double site_elevation,
                    double current_time, double *site_x, double *site_y,
                    double *site_z, double *site_vx, double *site_vy,
                    mat3x3 site_matrix)
{
   static double g1, g2; /* Flattening correction. */
   static double cos_lat, sin_lat;
   static double old_site_lat = -100000; /* Avoid unnecessary recompute. */
   static double old_site_elevation = -100000;
   double lat;
   double site_ra; /* Right ascension of site. */
   double cos_ra, sin_ra;

   if ((site_lat != old_site_lat) || (site_elevation != old_site_elevation)) {
      old_site_lat = site_lat;
      old_site_elevation = site_elevation;
      lat = atan(1 / (1 - SQR(EarthFlat)) * tan(site_lat));

      cos_lat = cos(lat);
      sin_lat = sin(lat);

      g1 = EarthRadius /
           (sqrt(1 - (2 * EarthFlat - SQR(EarthFlat)) * SQR(sin_lat)));
      g2 = g1 * SQR(1 - EarthFlat);
      g1 += site_elevation;
      g2 += site_elevation;
   }

   site_ra = PI2 * ((current_time - SidDay) * SiderealSolar + SidReference) -
             site_long;
   cos_ra = cos(site_ra);
   sin_ra = sin(site_ra);

   *site_x = g1 * cos_lat * cos_ra;
   *site_y = g1 * cos_lat * sin_ra;
   *site_z = g2 * sin_lat;
   *site_vx = -SidRate * (*site_y);
   *site_vy = SidRate * (*site_x);

   site_matrix[0][0] = sin_lat * cos_ra;
   site_matrix[0][1] = sin_lat * sin_ra;
   site_matrix[0][2] = -cos_lat;
   site_matrix[1][0] = -sin_ra;
   site_matrix[1][1] = cos_ra;
   site_matrix[1][2] = 0.0;
   site_matrix[2][0] = cos_ra * cos_lat;
   site_matrix[2][1] = sin_ra * cos_lat;
   site_matrix[2][2] = sin_lat;
}

/**
 * GetRange - compute range and range rate between site and satellite
 * @site_x: site x coordinate (km)
 * @site_y: site y coordinate (km)
 * @site_z: site z coordinate (km)
 * @site_vx: site x velocity (km/s)
 * @site_vy: site y velocity (km/s)
 * @sat_x: satellite x coordinate (km)
 * @sat_y: satellite y coordinate (km)
 * @sat_z: satellite z coordinate (km)
 * @sat_vx: satellite x velocity (km/s)
 * @sat_vy: satellite y velocity (km/s)
 * @sat_vz: satellite z velocity (km/s)
 * @range: output range (km)
 * @range_rate: output range rate (km/s)
 */
void GetRange(double site_x, double site_y, double site_z, double site_vx,
              double site_vy, double sat_x, double sat_y, double sat_z,
              double sat_vx, double sat_vy, double sat_vz, double *range,
              double *range_rate)
{
   double dx, dy, dz;

   dx = sat_x - site_x;
   dy = sat_y - site_y;
   dz = sat_z - site_z;

   *range = sqrt(SQR(dx) + SQR(dy) + SQR(dz));

   *range_rate =
      ((sat_vx - site_vx) * dx + (sat_vy - site_vy) * dy + sat_vz * dz) /
      (*range);
}

/**
 * GetTopocentric - convert geocentric coordinates to topocentric
 * @sat_x: satellite x coordinate (km)
 * @sat_y: satellite y coordinate (km)
 * @sat_z: satellite z coordinate (km)
 * @site_x: site x coordinate (km)
 * @site_y: site y coordinate (km)
 * @site_z: site z coordinate (km)
 * @site_matrix: transformation matrix produced by GetSitPosition()
 * @x: output topocentric x
 * @y: output topocentric y
 * @z: output topocentric z
 */
static void GetTopocentric(double sat_x, double sat_y, double sat_z,
                           double site_x, double site_y, double site_z,
                           mat3x3 site_matrix, double *x, double *y,
                           double *z)
{
   sat_x -= site_x;
   sat_y -= site_y;
   sat_z -= site_z;

   *x = site_matrix[0][0] * sat_x + site_matrix[0][1] * sat_y +
        site_matrix[0][2] * sat_z;
   *y = site_matrix[1][0] * sat_x + site_matrix[1][1] * sat_y +
        site_matrix[1][2] * sat_z;
   *z = site_matrix[2][0] * sat_x + site_matrix[2][1] * sat_y +
        site_matrix[2][2] * sat_z;
}

/**
 * GetBearings - compute azimuth/elevation bearings for a satellite
 * @sat_x: satellite x coordinate (km)
 * @sat_y: satellite y coordinate (km)
 * @sat_z: satellite z coordinate (km)
 * @site_x: site x coordinate (km)
 * @site_y: site y coordinate (km)
 * @site_z: site z coordinate (km)
 * @site_matrix: transformation matrix from GetSitPosition()
 * @azimuth: output azimuth (radians)
 * @elevation: output elevation (radians)
 */
void GetBearings(double sat_x, double sat_y, double sat_z, double site_x,
                 double site_y, double site_z, mat3x3 site_matrix,
                 double *azimuth, double *elevation)
{
   double x, y, z;

   GetTopocentric(sat_x, sat_y, sat_z, site_x, site_y, site_z, site_matrix,
                 &x, &y, &z);

   *elevation = atan(z / sqrt(SQR(x) + SQR(y)));
   *azimuth = PI - atan2(y, x);

   if (*azimuth < 0)
      *azimuth += PI;
}

/**
 * Eclipsed - determine whether a satellite is eclipsed by Earth
 * @sat_x: satellite x coordinate (km)
 * @sat_y: satellite y coordinate (km)
 * @sat_z: satellite z coordinate (km)
 * @sat_radius: satellite radius from geocenter (km)
 * @current_time: current time (days)
 *
 * Return: 1 if eclipsed, 0 otherwise.
 */
int Eclipsed(double sat_x, double sat_y, double sat_z, double sat_radius,
             double current_time)
{
   double mean_anomaly, true_anomaly;
   double sun_x, sun_y, sun_z, sun_rad;
   double vx, vy, vz;
   double cos_theta;

   mean_anomaly = SunMeanAnomaly +
                  (current_time - SunEpochTime) * SunMeanMotion * PI2;
   true_anomaly = Kepler(mean_anomaly, SunEccentricity);

   GetSatPosition(SunEpochTime, SunRAAN, SunArgPerigee, SunSemiMajorAxis,
                 SunInclination, SunEccentricity, 0.0, 0.0, current_time,
                 true_anomaly, &sun_x, &sun_y, &sun_z, &sun_rad, &vx, &vy,
                 &vz);

   cos_theta =
      (sun_x * sat_x + sun_y * sat_y + sun_z * sat_z) / (sun_rad * sat_radius) *
         CosPenumbra +
      (sat_radius / EarthRadius) * SinPenumbra;

   if (cos_theta < 0) {
      if (cos_theta <
          -sqrt(SQR(sat_radius) - SQR(EarthRadius)) / sat_radius * CosPenumbra +
             (sat_radius / EarthRadius) * SinPenumbra)
         return 1;
   }

   return 0;
}

/**
 * InitOrbitRoutines - initialize Sun ephemeris and sidereal reference
 * @epoch_day: epoch day number (days)
 *
 * Formulas are from "Explanatory Supplement to the Astronomical Ephemeris".
 */
void InitOrbitRoutines(double epoch_day)
{
   double t, t2, t3, omega;
   int n;
   double sun_true_anomaly, sun_distance;

   t = (floor(epoch_day) - 0.5) / 36525;
   t2 = t * t;
   t3 = t2 * t;

   SidDay = floor(epoch_day);

   SidReference = (6.6460656 + 2400.051262 * t + 0.00002581 * t2) / 24;
   SidReference -= floor(SidReference);

   /* Omega is used to correct for nutation and abberation. */
   omega = (259.18 - 1934.142 * t) * RadiansPerDegree;
   n = (int)(omega / PI2);
   omega -= n * PI2;

   SunEpochTime = epoch_day;
   SunRAAN = 0;

   SunInclination =
      (23.452294 - 0.0130125 * t - 0.00000164 * t2 + 0.000000503 * t3 +
       0.00256 * cos(omega)) *
      RadiansPerDegree;
   SunEccentricity = (0.01675104 - 0.00004180 * t - 0.000000126 * t2);
   SunArgPerigee =
      (281.220833 + 1.719175 * t + 0.0004527 * t2 + 0.0000033 * t3) *
      RadiansPerDegree;
   SunMeanAnomaly =
      (358.475845 + 35999.04975 * t - 0.00015 * t2 - 0.00000333333 * t3) *
      RadiansPerDegree;
   n = (int)(SunMeanAnomaly / PI2);
   SunMeanAnomaly -= n * PI2;

   SunMeanMotion = 1 / (365.24219879 - 0.00000614 * t);

   sun_true_anomaly = Kepler(SunMeanAnomaly, SunEccentricity);
   sun_distance = SunSemiMajorAxis * (1 - SQR(SunEccentricity)) /
                  (1 + SunEccentricity * cos(sun_true_anomaly));

   SinPenumbra = (SunRadius - EarthRadius) / sun_distance;
   CosPenumbra = sqrt(1 - SQR(SinPenumbra));
}

/**
 * SPrintTime - format a time value
 * @str: output string (must be large enough for HHMM:SS)
 * @time: time in days
 */
static void SPrintTime(char *str, double time)
{
   int day, hours, minutes, seconds;

   day = (int)time;
   time -= day;
   if (time < 0)
      time += 1.0; /* Correct for truncation problems with negatives. */

   hours = (int)(time * 24);
   time -= hours / 24.0;

   minutes = (int)(time * MinutesPerDay);
   time -= minutes / MinutesPerDay;

   seconds = (int)(time * SecondsPerDay);

   snprintf(str, 16, "%02d%02d:%02d", hours, minutes, seconds);
}

/**
 * PrintTime - print time as HHMM:SS
 * @out_file: output stream
 * @time: time in days
 */
void PrintTime(FILE *out_file, double time)
{
   char str[100];

   SPrintTime(str, time);
   fprintf(out_file, "%s", str);
}

/* Date calculation routines.
 * Robert Berger @ Carnegie Mellon
 *
 * January 1, 1900 is day 0.
 * Valid from 1900 through 2099.
 */
static const char *const MonthNames[] = {"Jan", "Feb", "Mar", "Apr", "May",
                                        "Jun", "Jul", "Aug", "Sep", "Oct",
                                        "Nov", "Dec"};

static const int MonthDays[] = {0, 31, 59, 90, 120, 151, 181,
                                212, 243, 273, 304, 334};

static const char *const DayNames[] = {"Sunday",   "Monday", "Tuesday",
                                      "Wednesday", "Thursday",
                                      "Friday",   "Saturday"};

/**
 * GetDayNum - get day number for a date
 * @year: year (2-digit values are interpreted as in the original code)
 * @month: month 1..12
 * @day: day-of-month (1..31)
 *
 * January 1 of the reference year is day 0.  Day number may be negative.
 *
 * Return: day number.
 */
long GetDayNum(int year, int month, int day)
{
   long result;

   /* Heuristic to allow 4 or 2 digit year specifications. */
   if (year < 50)
      year += 2000;
   else if (year < 100)
      year += 1900;

   result = ((((long)year - 1901) * 1461) >> 2) + MonthDays[month - 1] +
            day + 365;
   if (year % 4 == 0 && month > 2)
      result++;

   return result;
}

/**
 * GetDate - convert a day number to (year, month, day)
 * @day_num: day number
 * @yearp: output year
 * @monthp: output month (1..12)
 * @dayp: output day-of-month (1..31)
 */
static void GetDate(long day_num, int *yearp, int *monthp, int *dayp)
{
   int m;
   int l;
   long y;

   y = 4 * day_num;
   y /= 1461;

   day_num = day_num - 365 - (((y - 1) * 1461) >> 2);

   l = 0;
   if (y % 4 == 0 && day_num > MonthDays[2])
      l = 1;

   m = 1;
   while (day_num > MonthDays[m] + l)
      m++;

   day_num -= MonthDays[m - 1];
   if (m > 2)
      day_num -= l;

   *yearp = (int)(y + 1900);
   *monthp = m;
   *dayp = (int)day_num;
}

/**
 * GetDayOfWeek - get day-of-week for a day number
 * @day_num: day number
 *
 * Return: day-of-week (Sunday = 0).
 */
int GetDayOfWeek(long day_num)
{
   return (int)(day_num % 7);
}

/**
 * SPrintDate - format a date string
 * @str: output buffer
 * @day_num: day number
 */
static void SPrintDate(char *str, long day_num)
{
   int month, day, year;

   GetDate(day_num, &year, &month, &day);
   snprintf(str, 64, "%d %s %d", day, MonthNames[month - 1], year);
}

/**
 * SPrintDayOfWeek - format a day-of-week string
 * @str: output buffer
 * @day_num: day number
 */
static void SPrintDayOfWeek(char *str, long day_num)
{
   strcpy(str, DayNames[day_num % 7]);
}

/**
 * PrintDate - print date string
 * @out_file: output stream
 * @day_num: day number
 */
void PrintDate(FILE *out_file, long day_num)
{
   char str[100];

   SPrintDate(str, day_num);
   fprintf(out_file, "%s", str);
}

/**
 * PrintDayOfWeek - print day-of-week
 * @out_file: output stream
 * @day_num: day number
 */
void PrintDayOfWeek(FILE *out_file, long day_num)
{
   char str[100];

   SPrintDayOfWeek(str, day_num);
   fprintf(out_file, "%s", str);
}

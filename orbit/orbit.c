/* Copyright (c) 1986,1987,1988,1989,1990 Robert W. Berger N3EMO
   May be freely distributed, provided this notice remains intact. */

/* Change Log
   4/2/1990   v3.9 Misc bug fixes. Changed ElementSet and
         EpochRev to unsigned longs in nasa.c. Allow
         satellite names with spaces in mode.dat.

   3/15/1990   v3.8 Stop assigning single character abbreviations
         after the 62nd satellite.

   3/7/1990   v3.7 Make Phase III style phase (0-255) the default.
         Ignore case in satellite names.

   12/19/1989   v3.6 Use more direct calculations for dates.
         Calculate a new sidereal time reference for each run.

   12/8/1988   v3.5 Allow multiple overlapping modes in "mode.dat".

   6/28/1988   v3.4 Cleaned up Eclipse code. Fixed leap year handling
         for centesimal years. Added a heuristic to GetDay to
         allow 2 or 4 digit year specifications.

   1/25/1988   v3.2 Rewrote orbitr.c to improve modularity,
         efficiency, and accuracy. Adopted geocentric
         cartesian coordinates as the standard representation
         for position and velocity. Added direct calculation
          of range-rate for better doppler predections.

   12/1/1988   v3.1 Allow spaces in satellite names. Provide
         single character aliases for 62 satellites
         (up from 26).

   4/7/87      v3.0 Added eclipses.

   4/1/87      v2.4 Added "Flip" option in site file for
         0-180 elevation support.

   3/24/87      v2.3 Adapted for new kepler.dat format.
         Allow beacon frequencies in mode.dat.
         Use decay rate for drag compensation.

   5/10/86      v2.2 Added single character aliases for satellite
         names.

   4/30/86      v2.1 Print blank line if satellite dips below
         horizon and reappears during same orbit and day

   4/29/86      v2.0 Changed GetSatelliteParams() to use AMSAT's
         "kepler.dat" file. Moved schedule to "mode.dat" file.

   4/22/86      v1.3 Inserted N8FJB's suggestions for variable naming
         which maintain 8 character uniqueness.
         Also removed "include" file orbitr.h, which had two
         definitions of external functions defined in orbit.c
             -K3MC

   4/1/86      v1.2 Corrected a scanf conversion to %d for an int
         type.    -K3MC

   3/19/86      v1.1 Changed GetSatelliteParams to not pass NULL
         to sscanf.
                           */

#define DRAG 1

#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef double mat3x3[3][3];

/* External routines provided by orbitr.c. */
double Kepler(double mean_anomaly, double eccentricity);
int Eclipsed(double sat_x, double sat_y, double sat_z, double sat_radius,
             double current_time);
long GetDayNum(int year, int month, int day);
void GetBearings(double sat_x, double sat_y, double sat_z,
                 double site_x, double site_y, double site_z,
                 mat3x3 site_matrix, double *azimuth, double *elevation);
void GetPrecession(double semi_major_axis, double eccentricity,
                   double inclination, double *raan_precession,
                   double *perigee_precession);
void GetRange(double site_x, double site_y, double site_z,
              double site_vx, double site_vy,
              double sat_x, double sat_y, double sat_z,
              double sat_vx, double sat_vy, double sat_vz,
              double *range, double *range_rate);
void GetSatPosition(double epoch_time, double epoch_raan,
                    double epoch_arg_perigee, double semi_major_axis,
                    double inclination, double eccentricity,
                    double raan_precession, double perigee_precession,
                    double time, double true_anomaly,
                    double *x, double *y, double *z, double *radius,
                    double *vx, double *vy, double *vz);
void GetSitPosition(double site_lat, double site_long, double site_elevation,
                    double current_time, double *site_x, double *site_y,
                    double *site_z, double *site_vx, double *site_vy,
                    mat3x3 site_matrix);
void GetSubSatPoint(double sat_x, double sat_y, double sat_z, double time,
                    double *latitude, double *longitude, double *height);
void InitOrbitRoutines(double epoch_day);
void PrintDate(FILE *out_file, long day_num);
void PrintDayOfWeek(FILE *out_file, long day_num);
void PrintTime(FILE *out_file, double time);

#ifndef PI
#define PI 3.1415926535897932384626433832795028841971693993751058
#endif

#ifdef PI2
#undef PI2
#endif

#define PI2 (PI * 2)

#define MinutesPerDay (24 * 60.0)
#define SecondsPerDay (60 * MinutesPerDay)
#define HalfSecond (0.5 / SecondsPerDay)
#define EarthRadius 6378.16               /* Kilometers */
#define C 2.99792458e5                    /* Kilometers/Second */
#define DegreesPerRadian (180 / PI)
#define RadiansPerDegree (PI / 180)
#define ABS(x) ((x) < 0 ? (-(x)) : (x))
#define SQR(x) ((x) * (x))

#define MaxModes 10

typedef struct {
   int MinPhase;
   int MaxPhase;
   char ModeStr[20];
} ModeRec;

static const char VersionStr[] = "N3EMO Orbit Simulator  v3.9";

/* Keplerian elements and misc. data for the satellite. */
static double EpochDay;            /* time of epoch */
static double EpochMeanAnomaly;    /* Mean Anomaly at epoch */
static long EpochOrbitNum;         /* Integer orbit # of epoch */
static double EpochRAAN;           /* RAAN at epoch */
static double EpochMeanMotion;     /* Revolutions/day */
static double OrbitalDecay;        /* Revolutions/day^2 */
static double EpochArgPerigee;     /* argument of perigee at epoch */
static double Eccentricity;
static double Inclination;
static char SatName[100];
static int ElementSet;
static double BeaconFreq;          /* MHz, used for doppler calc */
static double MaxPhase;            /* Phase units in 1 orbit */
static double PerigeePhase;
static int NumModes;
static ModeRec Modes[MaxModes];
static int PrintApogee;
static int PrintEclipses;
static int Flip;

/* Simulation parameters. */
static double StartTime;
static double EndTime;
static double StepTime;            /* In days, 1 = New Year of reference year */

/* Site parameters. */
static char SiteName[100];
static double SiteLat;
static double SiteLong;
static double SiteAltitude;
static double SiteMinElev;

static int ListSatellites(void);
static void MatchStr(FILE *in_file, const char *file_name, const char *target);
static int LetterNum(char c);
static int cstrncmp(const char *str1, const char *str2, size_t len);
static void GetSatelliteParams(void);
static void GetSiteParams(void);
static void GetSimulationParams(void);
static void PrintMode(FILE *out_file, int phase);
static bool ReadLine(char *buf, size_t buf_len);
static void DiscardStdinLine(void);

/**
 * lc_char - ASCII-ish lowercase helper for case-insensitive compares
 * @c: input character
 *
 * This wraps ctype(3) to avoid undefined behavior on negative signed-char
 * values by casting through unsigned char.
 *
 * Return: lowercase form of @c if @c is an uppercase letter, otherwise @c.
 */
static inline int lc_char(int c)
{
   unsigned char uc = (unsigned char)c;

   if (isupper(uc))
      return tolower(uc);
   return (int)uc;
}

/**
 * ReadLine - read a line from stdin
 * @buf: destination buffer
 * @buf_len: size of @buf
 *
 * Reads a single line from stdin into @buf, strips trailing CR/LF, and
 * guarantees NUL-termination.
 *
 * Return: true on success, false on EOF/error.
 */
static bool ReadLine(char *buf, size_t buf_len)
{
   size_t n;

   if (buf_len == 0)
      return false;

   if (fgets(buf, (int)buf_len, stdin) == NULL)
      return false;

   n = strcspn(buf, "\r\n");
   buf[n] = '\0';
   return true;
}

/**
 * DiscardStdinLine - discard remainder of the current stdin line
 *
 * When mixing scanf() and fgets(), this is used to consume the trailing
 * newline left behind by scanf().
 */
static void DiscardStdinLine(void)
{
   int ch;

   do {
      ch = getchar();
   } while (ch != '\n' && ch != EOF);
}

/**
 * ListSatellites - list the satellites in kepler.dat
 *
 * Prints a catalog of "Satellite: ..." entries found in kepler.dat.
 *
 * Return: number of satellites found.
 */
static int ListSatellites(void)
{
   char str[100];
   FILE *in_file;
   char satchar;
   int num_satellites;

   printf("Available satellites:\n");

   in_file = fopen("kepler.dat", "r");
   if (in_file == NULL) {
      printf("\"kepler.dat\" not found\n");
      exit(EXIT_FAILURE);
   }

   satchar = 'a';
   num_satellites = 0;
   while (fgets(str, (int)sizeof(str), in_file) != NULL) {
      if (strncmp(str, "Satellite: ", 11) != 0)
         continue;

      printf("\t");
      if (satchar != 0)
         printf("%c) ", satchar);
      else
         printf("   ");

      printf("%s", &str[11]);

      if (satchar) {
         if (satchar == 'z')
            satchar = 'A';
         else if (satchar == 'Z')
            satchar = '0';
         else if (satchar == '9')
            satchar = 0;
         else
            satchar++;
      }

      num_satellites++;
   }

   fclose(in_file);
   return num_satellites;
}

/**
 * MatchStr - match and skip over a string in the input file
 * @in_file: input stream
 * @file_name: filename used in diagnostics
 * @target: string expected at the current file position
 *
 * Reads strlen(@target) bytes from @in_file and compares them against
 * @target.  Exits on mismatch.
 */
static void MatchStr(FILE *in_file, const char *file_name, const char *target)
{
   char str[100];
   size_t target_len = strlen(target);
   int nread;

   if (target_len + 1 > sizeof(str)) {
      printf("internal error: target too long in MatchStr()\n");
      exit(EXIT_FAILURE);
   }

   nread = (int)target_len + 1;
   if (fgets(str, nread, in_file) == NULL) {
      printf("%s: unexpected EOF while expecting \"%s\"\n", file_name,
             target);
      exit(EXIT_FAILURE);
   }

   if (strcmp(target, str) != 0) {
      printf("%s: found \"%s\" while expecting \"%s\"\n", file_name, str,
             target);
      exit(EXIT_FAILURE);
   }
}

/**
 * LetterNum - map a single-character satellite selector to a 1-based index
 * @c: selector character (a-z, A-Z, 0-9)
 *
 * Return: 1..62 on success, or 0 for invalid input.
 */
static int LetterNum(char c)
{
   if (c >= 'a' && c <= 'z')
      return (c - 'a') + 1;
   if (c >= 'A' && c <= 'Z')
      return (c - 'A') + 27;
   if (c >= '0' && c <= '9')
      return (c - '0') + 53;
   return 0;
}

/**
 * cstrncmp - case-insensitive strncmp
 * @str1: first string
 * @str2: second string
 * @len: number of characters to compare
 *
 * Return: 0 if equal (case-insensitive) for @len characters, 1 otherwise.
 */
static int cstrncmp(const char *str1, const char *str2, size_t len)
{
   size_t i;

   for (i = 0; i < len; i++) {
      if (lc_char((unsigned char)str1[i]) != lc_char((unsigned char)str2[i]))
         return 1;
   }

   return 0;
}

/**
 * GetSatelliteParams - prompt for satellite and load its elements
 *
 * Reads satellite selection from stdin, looks up the corresponding element
 * set in kepler.dat, and (optionally) loads mode information from mode.dat.
 */
static void GetSatelliteParams(void)
{
   FILE *in_file;
   char str[100];
   int epoch_year;
   int found;
   int i;
   int num_satellites;
   char satchar;
   int sat_index;

   num_satellites = ListSatellites();

   found = 0;
   while (!found) {
      printf("Letter or satellite name :");
      if (!ReadLine(SatName, sizeof(SatName)))
         exit(EXIT_FAILURE);

      in_file = fopen("kepler.dat", "r");
      if (in_file == NULL) {
         printf("kepler.dat not found\n");
         exit(EXIT_FAILURE);
      }

      if (strlen(SatName) == 1) {
         /* Use single character label. */
         satchar = SatName[0];
         sat_index = LetterNum(satchar);
         if (sat_index == 0 || sat_index > num_satellites) {
            printf("'%c' is out of range\n", satchar);
            fclose(in_file);
            continue;
         }

         for (i = 1; i <= sat_index; i++) {
            /* Find line beginning with "Satellite: ". */
            do {
               if (fgets(str, (int)sizeof(str), in_file) == NULL)
                  break;
            } while (strncmp(str, "Satellite: ", 11) != 0);
         }

         found = 1;
         {
            const char *start = &str[11];
            size_t n = strcspn(start, "\r\n");

            if (n >= sizeof(SatName))
               n = sizeof(SatName) - 1;
            memcpy(SatName, start, n);
            SatName[n] = '\0';
         }
      } else {
         /* Use satellite name. */
         while (!found) {
            if (fgets(str, (int)sizeof(str), in_file) == NULL)
               break; /* EOF */

            if (strncmp(str, "Satellite: ", 11) == 0) {
               if (cstrncmp(SatName, &str[11], strlen(SatName)) == 0)
                  found = 1;
            }
         }

         if (!found) {
            printf("Satellite %s not found\n", SatName);
            fclose(in_file);
         }
      }
   }

   BeaconFreq = 146.0; /* Default value (MHz). */

   /* Skip line following the Satellite: header. */
   if (fgets(str, (int)sizeof(str), in_file) == NULL)
      exit(EXIT_FAILURE);

   MatchStr(in_file, "kepler.dat", "Epoch time:");
   if (fgets(str, (int)sizeof(str), in_file) == NULL)
      exit(EXIT_FAILURE);
   if (sscanf(str, "%lf", &EpochDay) != 1)
      exit(EXIT_FAILURE);

   epoch_year = (int)(EpochDay / 1000.0);
   EpochDay -= epoch_year * 1000.0;
   EpochDay += (double)GetDayNum(epoch_year, 1, 0);
   if (fgets(str, (int)sizeof(str), in_file) == NULL)
      exit(EXIT_FAILURE);

   {
      long element_set;

      if (sscanf(str, "Element set: %ld", &element_set) != 1) {
         /* Old style kepler.dat */
         MatchStr(in_file, "kepler.dat", "Element set:");
         if (fgets(str, (int)sizeof(str), in_file) == NULL)
            exit(EXIT_FAILURE);
         if (sscanf(str, "%ld", &element_set) != 1)
            exit(EXIT_FAILURE);
      }
      ElementSet = (int)element_set;
   }

   MatchStr(in_file, "kepler.dat", "Inclination:");
   if (fgets(str, (int)sizeof(str), in_file) == NULL)
      exit(EXIT_FAILURE);
   if (sscanf(str, "%lf", &Inclination) != 1)
      exit(EXIT_FAILURE);
   Inclination *= RadiansPerDegree;

   MatchStr(in_file, "kepler.dat", "RA of node:");
   if (fgets(str, (int)sizeof(str), in_file) == NULL)
      exit(EXIT_FAILURE);
   if (sscanf(str, "%lf", &EpochRAAN) != 1)
      exit(EXIT_FAILURE);
   EpochRAAN *= RadiansPerDegree;

   MatchStr(in_file, "kepler.dat", "Eccentricity:");
   if (fgets(str, (int)sizeof(str), in_file) == NULL)
      exit(EXIT_FAILURE);
   if (sscanf(str, "%lf", &Eccentricity) != 1)
      exit(EXIT_FAILURE);

   MatchStr(in_file, "kepler.dat", "Arg of perigee:");
   if (fgets(str, (int)sizeof(str), in_file) == NULL)
      exit(EXIT_FAILURE);
   if (sscanf(str, "%lf", &EpochArgPerigee) != 1)
      exit(EXIT_FAILURE);
   EpochArgPerigee *= RadiansPerDegree;

   MatchStr(in_file, "kepler.dat", "Mean anomaly:");
   if (fgets(str, (int)sizeof(str), in_file) == NULL)
      exit(EXIT_FAILURE);
   if (sscanf(str, "%lf", &EpochMeanAnomaly) != 1)
      exit(EXIT_FAILURE);
   EpochMeanAnomaly *= RadiansPerDegree;

   MatchStr(in_file, "kepler.dat", "Mean motion:");
   if (fgets(str, (int)sizeof(str), in_file) == NULL)
      exit(EXIT_FAILURE);
   if (sscanf(str, "%lf", &EpochMeanMotion) != 1)
      exit(EXIT_FAILURE);

   MatchStr(in_file, "kepler.dat", "Decay rate:");
   if (fgets(str, (int)sizeof(str), in_file) == NULL)
      exit(EXIT_FAILURE);
   if (sscanf(str, "%lf", &OrbitalDecay) != 1)
      exit(EXIT_FAILURE);

   MatchStr(in_file, "kepler.dat", "Epoch rev:");
   if (fgets(str, (int)sizeof(str), in_file) == NULL)
      exit(EXIT_FAILURE);
   if (sscanf(str, "%ld", &EpochOrbitNum) != 1)
      exit(EXIT_FAILURE);

   while (1) {
      if (fgets(str, (int)sizeof(str), in_file) == NULL)
         break; /* EOF */
      if (strlen(str) <= 2)
         break; /* Blank line */
      (void)sscanf(str, "Beacon: %lf", &BeaconFreq);
   }

   PrintApogee = (Eccentricity >= 0.3);

   PerigeePhase = 0;
   MaxPhase = 256;
   NumModes = 0;

   fclose(in_file);

   in_file = fopen("mode.dat", "r");
   if (in_file == NULL)
      return;

   found = 0;
   while (!found) {
      if (fgets(str, (int)sizeof(str), in_file) == NULL)
         break; /* EOF */
      if (strncmp(str, "Satellite: ", 11) == 0) {
         if (cstrncmp(SatName, &str[11], strlen(SatName)) == 0)
            found = 1;
      }
   }

   if (found) {
      while (1) {
         if (fgets(str, (int)sizeof(str), in_file) == NULL)
            break; /* EOF */
         if (strlen(str) <= 2)
            break; /* Blank line */

         (void)sscanf(str, "Beacon: %lf", &BeaconFreq);
         (void)sscanf(str, "Perigee phase: %lf", &PerigeePhase);
         (void)sscanf(str, "Max phase: %lf", &MaxPhase);

         if (NumModes < MaxModes) {
            if (sscanf(str, "Mode: %19s from %d to %d",
                       Modes[NumModes].ModeStr,
                       &Modes[NumModes].MinPhase,
                       &Modes[NumModes].MaxPhase) == 3) {
               NumModes++;
            }
         }
      }
   }

   fclose(in_file);
}

/**
 * GetSiteParams - prompt for and load site parameters
 */
static void GetSiteParams(void)
{
   FILE *in_file;
   char name_in[100];
   char name[100];
   char str[100];

   printf("Site name :");
   if (!ReadLine(name_in, sizeof(name_in)))
      exit(EXIT_FAILURE);

   if (snprintf(name, sizeof(name), "%s.sit", name_in) >= (int)sizeof(name)) {
      printf("site name too long\n");
      exit(EXIT_FAILURE);
   }

   in_file = fopen(name, "r");
   if (in_file == NULL) {
      printf("%s not found\n", name);
      exit(EXIT_FAILURE);
   }

   if (fgets(SiteName, (int)sizeof(SiteName), in_file) == NULL)
      exit(EXIT_FAILURE);

   if (fgets(str, (int)sizeof(str), in_file) == NULL)
      exit(EXIT_FAILURE);
   (void)sscanf(str, "%lf", &SiteLat);
   SiteLat *= RadiansPerDegree;

   if (fgets(str, (int)sizeof(str), in_file) == NULL)
      exit(EXIT_FAILURE);
   (void)sscanf(str, "%lf", &SiteLong);
   SiteLong *= RadiansPerDegree;

   if (fgets(str, (int)sizeof(str), in_file) == NULL)
      exit(EXIT_FAILURE);
   (void)sscanf(str, "%lf", &SiteAltitude);
   SiteAltitude /= 1000; /* convert to km */

   if (fgets(str, (int)sizeof(str), in_file) == NULL)
      exit(EXIT_FAILURE);
   (void)sscanf(str, "%lf", &SiteMinElev);
   SiteMinElev *= RadiansPerDegree;

   Flip = 0;
   PrintEclipses = 0;
   while (fgets(str, (int)sizeof(str), in_file) != NULL) {
      if (strncmp(str, "Flip", 4) == 0)
         Flip = 1;
      else if (strncmp(str, "Eclipse", 7) == 0)
         PrintEclipses = 1;
      else
         printf("\"%s\" unknown option: %s", name, str);
   }

   fclose(in_file);
}

/**
 * GetSimulationParams - prompt for simulation start/end and time step
 */
static void GetSimulationParams(void)
{
   double hour;
   double duration;
   int month;
   int day;
   int year;

again:
   printf("Start date (UTC) (Month Day Year) :");
   if (scanf("%d%d%d", &month, &day, &year) != 3)
      exit(EXIT_FAILURE);
   if (year < 2000)
      goto again;

   StartTime = (double)GetDayNum(year, month, day);
   printf("Starting Hour (UTC) :");
   if (scanf("%lf", &hour) != 1)
      exit(EXIT_FAILURE);
   StartTime += hour / 24;

   printf("Duration (Days) :");
   if (scanf("%lf", &duration) != 1)
      exit(EXIT_FAILURE);
   EndTime = StartTime + duration;

   printf("Time Step (Minutes) :");
   if (scanf("%lf", &StepTime) != 1)
      exit(EXIT_FAILURE);
   StepTime /= MinutesPerDay;

   DiscardStdinLine();
}

/**
 * PrintMode - emit mode(s) active at a given phase
 * @out_file: output stream
 * @phase: current satellite phase
 */
static void PrintMode(FILE *out_file, int phase)
{
   int cur_mode;

   for (cur_mode = 0; cur_mode < NumModes; cur_mode++) {
      if ((phase >= Modes[cur_mode].MinPhase && phase < Modes[cur_mode].MaxPhase)
          || ((Modes[cur_mode].MinPhase > Modes[cur_mode].MaxPhase)
              && (phase >= Modes[cur_mode].MinPhase
                  || phase < Modes[cur_mode].MaxPhase))) {
         fprintf(out_file, "%s ", Modes[cur_mode].ModeStr);
      }
   }
}

/**
 * main - entry point
 *
 * Return: 0 on success, non-zero on failure.
 */
int main(void)
{
   double reference_orbit;    /* Floating point orbit # at epoch */
   double current_time;       /* In days */
   double current_orbit;
   double average_motion;     /* Corrected for drag */
   double current_motion;
   double mean_anomaly;
   double true_anomaly;
   double semi_major_axis;
   double radius;             /* From geocenter */
   double sat_x, sat_y, sat_z;
   double sat_vx, sat_vy, sat_vz; /* Kilometers/second */
   double site_x, site_y, site_z;
   double site_vx, site_vy;
   mat3x3 site_matrix;
   double height;
   double raan_precession;
   double perigee_precession;
   double ssp_lat;
   double ssp_long;
   long orbit_num, prev_orbit_num;
   long day_num, prev_day_num;
   double azimuth, elevation, range;
   double range_rate, doppler;
   int phase;
   char filename[100];
   FILE *out_file;
   int did_apogee;
   double tmp_time;
   int prev_visible;

   printf("%s\n", VersionStr);

   GetSatelliteParams();
   GetSiteParams();
   GetSimulationParams();

   InitOrbitRoutines((StartTime + EndTime) / 2);

   printf("Output file (RETURN for TTY) :");
   if (!ReadLine(filename, sizeof(filename)))
      exit(EXIT_FAILURE);

   if (strlen(filename) > 0) {
      out_file = fopen(filename, "w");
      if (out_file == NULL) {
         printf("Can't write to %s\n", filename);
         exit(EXIT_FAILURE);
      }
   } else {
      out_file = stdout;
   }

   fprintf(out_file, "%s Element Set %d\n", SatName, ElementSet);
   fprintf(out_file, "%s\n", SiteName);
   fprintf(out_file, "Doppler calculated for freq = %lf MHz\n", BeaconFreq);

   semi_major_axis = 331.25 * exp(2 * log(MinutesPerDay / EpochMeanMotion) / 3);
   GetPrecession(semi_major_axis, Eccentricity, Inclination,
                 &raan_precession, &perigee_precession);

   reference_orbit = EpochMeanAnomaly / PI2 + EpochOrbitNum;

   prev_day_num = -10000;
   prev_orbit_num = -10000;
   prev_visible = 0;

   BeaconFreq *= 1E6; /* Convert to Hz. */

   did_apogee = 0;

   for (current_time = StartTime; current_time <= EndTime;
        current_time += StepTime) {
      average_motion = EpochMeanMotion + (current_time - EpochDay) * OrbitalDecay / 2;
      current_motion = EpochMeanMotion + (current_time - EpochDay) * OrbitalDecay;

      semi_major_axis = 331.25 * exp(2 * log(MinutesPerDay / current_motion) / 3);

      current_orbit = reference_orbit + (current_time - EpochDay) * average_motion;
      orbit_num = (long)current_orbit;

      mean_anomaly = (current_orbit - (double)orbit_num) * PI2;

      tmp_time = current_time;
      if (mean_anomaly < PI)
         did_apogee = 0;
      if (PrintApogee && !did_apogee && mean_anomaly > PI) {
         /* Calculate apogee. */
         tmp_time -= StepTime; /* Pick up later where we left off. */
         mean_anomaly = PI;
         current_time = EpochDay +
                        (orbit_num - reference_orbit + 0.5) / average_motion;
      }

      true_anomaly = Kepler(mean_anomaly, Eccentricity);

      GetSatPosition(EpochDay, EpochRAAN, EpochArgPerigee, semi_major_axis,
                     Inclination, Eccentricity,
                     raan_precession, perigee_precession,
                     current_time, true_anomaly,
                     &sat_x, &sat_y, &sat_z, &radius,
                     &sat_vx, &sat_vy, &sat_vz);

      GetSitPosition(SiteLat, SiteLong, SiteAltitude, current_time,
                     &site_x, &site_y, &site_z, &site_vx, &site_vy,
                     site_matrix);

      GetBearings(sat_x, sat_y, sat_z, site_x, site_y, site_z, site_matrix,
                  &azimuth, &elevation);

      if (elevation >= SiteMinElev && current_time >= StartTime) {
         day_num = (long)(current_time + HalfSecond);
         if (((double)day_num) > current_time + HalfSecond)
            day_num -= 1; /* Correct for truncation of negative values. */

         if (orbit_num == prev_orbit_num && day_num == prev_day_num && !prev_visible)
            fprintf(out_file, "\n"); /* Dipped out of sight; print blank. */

         if (orbit_num != prev_orbit_num || day_num != prev_day_num) {
            /* Print header. */
            PrintDayOfWeek(out_file, day_num);
            fprintf(out_file, " ");
            PrintDate(out_file, day_num);
            fprintf(out_file, "  ----Orbit # %ld-----\n", orbit_num);
            fprintf(out_file, " U.T.C.   Az  El  ");
            if (Flip)
               fprintf(out_file, " Az'  El' ");

            fprintf(out_file, "Doppler Range");
            fprintf(out_file, " Height  Lat  Long  Phase(%3.0lf)\n", MaxPhase);
         }

         prev_orbit_num = orbit_num;
         prev_day_num = day_num;
         PrintTime(out_file, current_time + HalfSecond);

         fprintf(out_file, "  %3.0lf %3.0lf", azimuth * DegreesPerRadian,
                 elevation * DegreesPerRadian);

         if (Flip) {
            azimuth += PI;
            if (azimuth >= PI2)
               azimuth -= PI2;
            elevation = PI - elevation;
            fprintf(out_file, "  %3.0lf  %3.0lf", azimuth * DegreesPerRadian,
                    elevation * DegreesPerRadian);
         }

         GetRange(site_x, site_y, site_z, site_vx, site_vy,
                  sat_x, sat_y, sat_z, sat_vx, sat_vy, sat_vz,
                  &range, &range_rate);
         doppler = -BeaconFreq * range_rate / C;
         fprintf(out_file, "  %6.0lf %6.0lf", doppler, range);

         GetSubSatPoint(sat_x, sat_y, sat_z, current_time,
                        &ssp_lat, &ssp_long, &height);
         fprintf(out_file, " %6.0lf  %3.0lf  %4.0lf",
                 height,
                 ssp_lat * DegreesPerRadian,
                 ssp_long * DegreesPerRadian);

         phase = (int)(mean_anomaly / PI2 * MaxPhase + PerigeePhase);
         while (phase < 0)
            phase += (int)MaxPhase;
         while (phase >= (int)MaxPhase)
            phase -= (int)MaxPhase;

         fprintf(out_file, " %4d  ", phase);
         PrintMode(out_file, phase);

         if (PrintApogee && (mean_anomaly == PI))
            fprintf(out_file, "    Apogee");

         if (PrintEclipses) {
            if (Eclipsed(sat_x, sat_y, sat_z, radius, current_time))
               fprintf(out_file, "  Eclipse");
         }

         fprintf(out_file, "\n");
         prev_visible = 1;
      } else {
         prev_visible = 0;
      }

      if (PrintApogee && (mean_anomaly == PI))
         did_apogee = 1;

      current_time = tmp_time;
   }

   if (out_file != stdout)
      fclose(out_file);
   return 0;
}

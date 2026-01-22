/* nasa2kep.c - convert file of NASA keplerians to AMSAT format
 *
 * 7/12/87  Robert W. Berger N3EMO
 *
 * You know, if "0<space>" preceeds satellite name, just take it off
 * since I've seen files both ways.
 *
 * Mon Jun 17 15:50:22 PDT 2013
 *
 * 11/30/88  v3.5  Incorporate Greg Troxel's fix for reading epoch times
 *                 with imbedded spaces
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {
   BUF_LEN = 100
};

/**
 * StripLeadingZeroSpace - drop an optional "0 " prefix from a name line
 * @sat_name: Satellite name line buffer (NUL-terminated, may include '\n')
 *
 * Some TLE dumps prefix the satellite name line with "0 ".  Historically this
 * converter accepted both formats; preserve that behavior by stripping the
 * prefix when present.
 */
static void
StripLeadingZeroSpace(char *sat_name)
{
   if (sat_name[0] == '0' && sat_name[1] == ' ')
      memmove(sat_name, sat_name + 2, strlen(sat_name + 2) + 1);
}

/**
 * ReadOrDie - read a single text line or terminate with an error
 * @fp: Input stream
 * @buf: Destination buffer
 * @buf_len: Size of @buf in bytes
 * @what: Short label for diagnostics (e.g., "TLE line 1")
 */
static void
ReadOrDie(FILE *fp, char *buf, size_t buf_len, const char *what)
{
   if (fgets(buf, buf_len, fp) == NULL) {
      printf("Unexpected EOF while reading %s\n", what);
      exit(-1);
   }
}

/**
 * main - convert "nasa.dat" (NASA keplerians / TLE) to "kepler.dat" (AMSAT)
 *
 * This is a small conversion utility with fixed input/output file names, as in
 * the original implementation.
 */
int
main(void)
{
   char sat_name[BUF_LEN];
   char line1[BUF_LEN];
   char line2[BUF_LEN];
   FILE *in_file;
   FILE *out_file;
   int line_num;
   int sat_num;
   unsigned long element_set;
   unsigned long epoch_rev;
   double epoch_day;
   double decay_rate;
   double inclination;
   double raan;
   double eccentricity;
   double arg_perigee;
   double mean_anomaly;
   double mean_motion;
   double epoch_year;

   in_file = fopen("nasa.dat", "r");
   if (in_file == NULL) {
      printf("\"nasa.dat\" not found\n");
      exit(-1);
   }

   out_file = fopen("kepler.dat", "w");
   if (out_file == NULL) {
      printf("Can't write \"kepler.dat\"\n");
      fclose(in_file);
      exit(-1);
   }

   printf("Input File \"nasa.dat\"    Output File \"kepler.dat\"\n\n");
   /* Just a bit of user doc when one is rusty by not using the program much */

   sat_name[0] = sat_name[1] = sat_name[2] = (char)0;

   while (fgets(sat_name, sizeof(sat_name), in_file) != NULL) {
      if (sat_name[0] == '0')
         StripLeadingZeroSpace(sat_name);

      /* Squash out the "0 " that sometimes appears for the SatName */
      printf("%s", sat_name);

      ReadOrDie(in_file, line1, sizeof(line1), "TLE line 1");
      ReadOrDie(in_file, line2, sizeof(line2), "TLE line 2");

      if (sscanf(line1, "%1d", &line_num) != 1 || line_num != 1) {
         printf("Line 1 not present for satellite %s", sat_name);
         fclose(in_file);
         fclose(out_file);
         exit(-1);
      }

      if (sscanf(line2, "%1d", &line_num) != 1 || line_num != 2) {
         printf("Line 2 not present for satellite %s", sat_name);
         fclose(in_file);
         fclose(out_file);
         exit(-1);
      }

      if (sscanf(line1,
         "%*2c%5d%*10c%2lf%12lf%10lf%*21c%5lu",
         &sat_num, &epoch_year, &epoch_day, &decay_rate, &element_set) != 5) {
         printf("Failed to parse TLE line 1 for satellite %s", sat_name);
         fclose(in_file);
         fclose(out_file);
         exit(-1);
      }

      epoch_day += epoch_year * 1000.0;
      element_set /= 10; /* strip off checksum */

      if (sscanf(line2,
         "%*8c%8lf%8lf%7lf%8lf%8lf%11lf%5lu",
         &inclination, &raan, &eccentricity, &arg_perigee, &mean_anomaly, &mean_motion, &epoch_rev) != 7) {
         printf("Failed to parse TLE line 2 for satellite %s", sat_name);
         fclose(in_file);
         fclose(out_file);
         exit(-1);
      }

      epoch_rev /= 10; /* strip off checksum */
      eccentricity *= 1E-7;

      fprintf(out_file, "Satellite: %s", sat_name);
      fprintf(out_file, "Catalog number: %d\n", sat_num);
      fprintf(out_file, "Epoch time: %lf\n", epoch_day);
      fprintf(out_file, "Element set: %lu\n", element_set);
      fprintf(out_file, "Inclination: %lf deg\n", inclination);
      fprintf(out_file, "RA of node: %lf deg\n", raan);
      fprintf(out_file, "Eccentricity: %lf\n", eccentricity);
      fprintf(out_file, "Arg of perigee: %lf deg\n", arg_perigee);
      fprintf(out_file, "Mean anomaly: %lf deg\n", mean_anomaly);
      fprintf(out_file, "Mean motion: %lf rev/day\n", mean_motion);
      fprintf(out_file, "Decay rate: %le rev/day^2\n", decay_rate);
      fprintf(out_file, "Epoch rev: %lu\n", epoch_rev);
      fprintf(out_file, "\n");
   }

   fclose(in_file);
   fclose(out_file);

   return 0;
}

using System;
using System.Collections.Generic;
using System.Text;

namespace GeoRSS
{
    public class GeoUtils
    {
        public String EastingValue = "";
        public String NorthingValue = "";
        public String LatitudeValue = "";
        public String LongtitudeValue = "";
        public double deg2rad = Math.PI / 180;
        public double rad2deg = 180.0 / Math.PI;
        public double pi = Math.PI;

        public void GPSConversion(String GPSLocation)
        {
            string lat = "";
            string lon = "";
            string[] geo = GPSLocation.Split(',');
            if (geo.Length == 2)
            {
                lat = geo[0];
                lon = geo[1];
                val_lldecimal(lat, lon);
            }
        }
        public void ne2ll(String easting, String northing, int grid, int locres)  // (east, north, grid)
        {
            // converts NGR easting and nothing to lat, lon.
            // input metres, output radians
            //String grid = 1;
            double a = 0;
            double b = 0;
            int e0 = 0;
            int n0 = 0;
            double f0 = 0;
            double e2 = 0;
            double lam0 = 0;
            double phi0 = 0;


            double north = double.Parse(northing);
            double east = double.Parse(easting);
            //var locres = 6;
            if (grid == 2)  // Irish
            {
                a = 6377340.189;       // OSGBI semi-major
                b = 6356034.447;        // OSGBI semi-minor
                e0 = 200000;           // easting of false origin
                n0 = 250000;          // northing of false origin
                f0 = 1.000035;     // OSGBI scale factor on central meridian
                e2 = 0.00667054015;  // OSGBI eccentricity squared
                lam0 = -0.13962634015954636615389526147909;  // OSGBI false east
                phi0 = 0.93375114981696632365417456114141;    // OSGBI false north
            }
            if (grid == 1)  // British
            {
                a = 6377563.396;       // OSI semi-major
                b = 6356256.91;        // OSI semi-minor
                e0 = 400000;           // easting of false origin
                n0 = -100000;          // northing of false origin
                f0 = 0.9996012717;     // OSI scale factor on central meridian
                e2 = 0.0066705397616;  // OSI eccentricity squared
                lam0 = -0.034906585039886591;  // OSI false east
                phi0 = 0.85521133347722145;    // OSI false north
            }
            if (grid == 3)   // Channel Is
            {
                a = 6378388.000;       // INT24 ED50 semi-major
                b = 6356911.946;       // INT24 ED50 semi-minor 
                e0 = 500000;           // easting of false origin
                n0 = 0;                // northing of false origin
                f0 = 0.9996;           // INT24 ED50 scale factor on central meridian
                e2 = 0.0067226700223333;  // INT24 ED50 eccentricity squared
                lam0 = -0.0523598775598;  // INT24 ED50 false east
                phi0 = 0 * deg2rad;    // INT24 ED50 false north
            }
            double af0 = a * f0;
            double bf0 = b * f0;
            double n = (af0 - bf0) / (af0 + bf0);
            double Et = east - e0;
            double phid = InitialLat(north, n0, af0, phi0, n, bf0);
            double nu = af0 / (sqrt(1 - (e2 * (sin(phid) * sin(phid)))));
            double rho = (nu * (1 - e2)) / (1 - (e2 * (sin(phid)) * (sin(phid))));
            double eta2 = (nu / rho) - 1;
            double tlat2 = tan(phid) * tan(phid);
            double tlat4 = pow(tan(phid), 4);
            double tlat6 = pow(tan(phid), 6);
            double clatm1 = pow(cos(phid), -1);
            double VII = tan(phid) / (2 * rho * nu);
            double VIII = (tan(phid) / (24 * rho * (nu * nu * nu))) * (5 + (3 * tlat2) + eta2 - (9 * eta2 * tlat2));
            double IX = ((tan(phid)) / (720 * rho * pow(nu, 5))) * (61 + (90 * tlat2) + (45 * tlat4));
            double phip = (phid - ((Et * Et) * VII) + (pow(Et, 4) * VIII) - (pow(Et, 6) * IX));
            double X = pow(cos(phid), -1) / nu;
            double XI = (clatm1 / (6 * (nu * nu * nu))) * ((nu / rho) + (2 * (tlat2)));
            double XII = (clatm1 / (120 * pow(nu, 5))) * (5 + (28 * tlat2) + (24 * tlat4));
            double XIIA = clatm1 / (5040 * pow(nu, 7)) * (61 + (662 * tlat2) + (1320 * tlat4) + (720 * tlat6));
            double lambdap = (lam0 + (Et * X) - ((Et * Et * Et) * XI) + (pow(Et, 5) * XII) - (pow(Et, 7) * XIIA));
            double latitude = 0;
            double longitude = 0;
            convert_to_wgs(grid, phip, lambdap, out latitude, out longitude);
            double lat = latitude * rad2deg;
            double lon = longitude * rad2deg;
            LatitudeValue = lat.ToString();
            LongtitudeValue = lon.ToString();
            //alert(lat);
            //alert(lon);
            todms(lat, lon);
            //var loc = calcloc(lon, lat);
            //if (locres == -1)                // safari bug -1 after reset.
            //locres = 0;
            //var locy = loc.substring(0, 6 + 2 * locres);
            //alert(locy);
        }
        public void todms(double lat, double lon)
        {
            int latbrg = 1;
            int lonbrg = 2;
            if (lat < 0)
                latbrg = 2;
            if (lon < 0)
                lonbrg = 1;
            double tlat = abs(lat);
            double tlon = abs(lon);
            double deglat = floor(tlat);
            double t = (tlat - deglat) * 60;
            double minlat2 = floor((t * 1000) + 0.0005);  //mins and hundredths ****
            minlat2 = minlat2 / 1000;
            double minlat = floor(t);
            double seclat = (t - minlat) * 60;
            seclat = seclat + 0.005;  //round to .005 sec
            seclat = floor(seclat * 100) / 100;  //  works in js 1.4
            double deglon = floor(tlon);
            t = (tlon - deglon) * 60;
            double minlon = floor(t);
            double minlon2 = floor((t * 1000) + 0.0005);  //mins and hundredths ****
            minlon2 = minlon2 / 1000;
            double seclon = (t - minlon) * 60;
            seclon = seclon + 0.005;
            seclon = floor(seclon * 100) / 100;  // js 1.4

            //MessageBox.Show(deglat.ToString());
            //MessageBox.Show(minlat.ToString());
            //MessageBox.Show(minlat2.ToString());
            //MessageBox.Show(seclat.ToString());
            //MessageBox.Show(latbrg.ToString());
            //MessageBox.Show(deglon.ToString());
            //MessageBox.Show(minlon.ToString());
            //MessageBox.Show(minlon2.ToString());
            //MessageBox.Show(seclon.ToString());
            //MessageBox.Show(lonbrg.ToString());
        }
        public double InitialLat(double north, double n0, double af0, double phi0, double n, double bf0)
        {
            double phi1 = ((north - n0) / af0) + phi0;
            double M = Marc(bf0, n, phi0, phi1);
            double phi2 = ((north - n0 - M) / af0) + phi1;
            double ind = 0;
            while ((abs(north - n0 - M) > 0.00001) && (ind < 20))  // max 20 iterations in case of error
            {
                ind = ind + 1;
                phi2 = ((north - n0 - M) / af0) + phi1;
                M = Marc(bf0, n, phi0, phi2);
                phi1 = phi2;
            }
            return (phi2);
        }
        public void val_ll(String latd, String latm, String lats, String latbs, String lond, String lonm, String lons, String lonbs, String locres)
        {
            int latb = 0;
            int lonb = 0;
            String grid = "British";

            if (latd == "")
            {
                //err.Text = "Latitude degrees must be entered";
                return;
            }

            if (latbs == "N")
                latb = 1;
            else if (latbs == "S")
                latb = 2;

            if (lonbs == "W")
                lonb = 1;
            else if (lonbs == "E")
                lonb = 2;

            if (abs(Number(latd)) >= 90)
            {
                //err.Text = "Degrees wrong";
                return;
            }
            if (Number(latm) >= 60)
            {
                //err.Text = "Minutes wrong";
                return;
            }
            if (Number(lats) >= 60)
            {
                //err.Text = "Seconds wrong";
                return;
            }
            if (abs(Number(lond)) >= 180)
            {
                //err.Text = "Degrees wrong";
                return;
            }
            if (Number(lonm) >= 60)
            {
                //err.Text = "Minutes wrong";
                return;
            }
            if (Number(lons) >= 60)
            {
                //err.Text = "Seconds wrong";
                return;
            }

            double lat = Number(latd);
            lat = lat + Number(latm) / 60;
            lat = lat + Number(lats) / 3600;
            if (latb == 2)  // S
                lat = lat * -1;
            double lon = Number(lond);
            lon = lon + Number(lonm) / 60;
            lon = lon + Number(lons) / 3600;
            if (lonb == 1)  // W
                lon = lon * -1;
            double latm2 = Number(latm) + Number(lats) / 60;
            latm2 = floor(latm2 * 1000) / 1000;
            double lonm2 = Number(lonm) + Number(lons) / 60;
            lonm2 = floor(lonm2 * 1000) / 1000;

            //double loc = calcloc(lon, lat);  //wgs84
            //int locy = loc.substring(0, 6 + 2 * locres);
            //int grid = choose_where(lat, lon);
            double phip = lat * deg2rad;  // deg to rad
            double lambdap = lon * deg2rad;
            //if (grid == "")
            //    return;
            double latitude = 0;
            double longitude = 0;
            convert_to_local(grid, phip, lambdap, out latitude, out longitude);
            phip = latitude;
            lambdap = longitude;

            ll2ne(phip, lambdap, grid);
        }
        public void val_lldecimal(String lats, String lons)
        {
            double lat = Number(lats);
            double lon = Number(lons);
            String grid = "British";
            double phip = lat * deg2rad;  // deg to rad
            double lambdap = lon * deg2rad;
            //if (grid == "")
            //    return;
            double latitude = 0;
            double longitude = 0;
            convert_to_local(grid, phip, lambdap, out latitude, out longitude);
            phip = latitude;
            lambdap = longitude;

            ll2ne(phip, lambdap, grid);
        }
        public void convert_to_local(String grid, double phip, double lambdap, out double latitude, out double longitude)
        {
            double WGS84_AXIS = 6378137;
            double WGS84_ECCENTRIC = 0.00669438037928458;
            double OSGB_AXIS = 6377563.396;
            double OSGB_ECCENTRIC = 0.0066705397616;
            double IRISH_AXIS = 6377340.189;
            double IRISH_ECCENTRIC = 0.00667054015;
            double INT24_AXIS = 6378388.000;
            double INT24_ECCENTRIC = 0.0067226700223333;
            double height = 10;  // dummy height

            double lat = 0;
            double lon = 0;

            if (grid == "British")
            {
                transform(phip, lambdap, WGS84_AXIS, WGS84_ECCENTRIC, height, OSGB_AXIS, OSGB_ECCENTRIC, -446.448, 125.157, -542.06, -0.1502, -0.247, -0.8421, 20.4894, out lat, out lon);
            }
            if (grid == "Irish")
            {
                transform(phip, lambdap, WGS84_AXIS, WGS84_ECCENTRIC, height, IRISH_AXIS, IRISH_ECCENTRIC, -482.53, 130.596, -564.557, 1.042, 0.214, 0.631, -8.15, out lat, out lon);
            }
            if (grid == "Channel Islands")
            {
                transform(phip, lambdap, WGS84_AXIS, WGS84_ECCENTRIC, height, INT24_AXIS, INT24_ECCENTRIC, 83.901, 98.127, 118.635, 0, 0, 0, 0, out lat, out lon);
            }
            latitude = lat;
            longitude = lon;
        }
        public void convert_to_wgs(int grid, double phip, double lambdap, out double latitude, out double longitude)
        {
            double WGS84_AXIS = 6378137;
            double WGS84_ECCENTRIC = 0.00669438037928458;
            double OSGB_AXIS = 6377563.396;
            double OSGB_ECCENTRIC = 0.0066705397616;
            double IRISH_AXIS = 6377340.189;
            double IRISH_ECCENTRIC = 0.00667054015;
            double INT24_AXIS = 6378388.000;
            double INT24_ECCENTRIC = 0.0067226700223333;
            double height = 10;  // dummy height

            double lat = 0;
            double lon = 0;

            if (grid == 1)
            {
                transform(phip, lambdap, OSGB_AXIS, OSGB_ECCENTRIC, height, WGS84_AXIS, WGS84_ECCENTRIC, 446.448, -125.157, 542.06, 0.1502, 0.247, 0.8421, 20.4894, out lat, out lon);
            }
            if (grid == 2)
            {
                transform(phip, lambdap, IRISH_AXIS, IRISH_ECCENTRIC, height, WGS84_AXIS, WGS84_ECCENTRIC, 482.53, -130.596, 564.557, -1.042, -0.214, -0.631, -8.15, out lat, out lon);
            }
            if (grid == 3)
            {
                transform(phip, lambdap, INT24_AXIS, INT24_ECCENTRIC, height, WGS84_AXIS, WGS84_ECCENTRIC, -83.901, -98.127, -118.635, 0, 0, 0, 0, out lat, out lon);
            }

            latitude = lat;
            longitude = lon;
        }
        public void transform(double lat, double lon, double a, double e, double h, double a2, double e2, double xp, double yp, double zp, double xr, double yr, double zr, double s, out double latitude, out double longitude)
        {
            // convert to cartesian; lat, lon are radians
            double sf = s * 0.000001;
            double v = a / (sqrt(1 - (e * (sin(lat) * sin(lat)))));
            double x = (v + h) * cos(lat) * cos(lon);
            double y = (v + h) * cos(lat) * sin(lon);
            double z = ((1 - e) * v + h) * sin(lat);

            double xrot = (xr / 3600) * deg2rad;
            double yrot = (yr / 3600) * deg2rad;
            double zrot = (zr / 3600) * deg2rad;

            double hx = x + (x * sf) - (y * zrot) + (z * yrot) + xp;
            double hy = (x * zrot) + y + (y * sf) - (z * xrot) + yp;
            double hz = (-1 * x * yrot) + (y * xrot) + z + (z * sf) + zp;

            // Convert back to lat, lon
            lon = atan(hy / hx);
            double p = sqrt((hx * hx) + (hy * hy));
            lat = atan(hz / (p * (1 - e2)));
            v = a2 / (sqrt(1 - e2 * (sin(lat) * sin(lat))));
            double errvalue = 1.0;
            double lat0 = 0;
            while (errvalue > 0.001)
            {
                lat0 = atan((hz + e2 * v * sin(lat)) / p);
                errvalue = abs(lat0 - lat);
                lat = lat0;
            }
            h = p / cos(lat) - v;
            latitude = lat;
            longitude = lon;
            //var geo = { latitude: lat, longitude: lon };
            //return(geo);
        }
        public void ll2ne(double lat, double lon, String grid)
        {
            double a = 0;
            double b = 0;
            int e0 = 0;
            int n0 = 0;
            double f0 = 0;
            double e2 = 0;
            double lam0 = 0;
            double phi0 = 0;

            // converts Lat/Lon OSGB to Easting and Northing: inputs radians.
            double phi = lat;   // latitude rad
            double lam = lon;   // longitude rad
            // grid = grid;       // British, Irish, Channel
            if (grid == "British")
            {
                a = 6377563.396;       // OSGB semi-major
                b = 6356256.91;        // OSGB semi-minor
                e0 = 400000;           // easting of false origin
                n0 = -100000;          // northing of false origin
                f0 = 0.9996012717;     // OSGB scale factor on central meridian
                e2 = 0.0066705397616;  // OSGB eccentricity squared
                lam0 = -0.034906585039886591;  // OSGB false east
                phi0 = 0.85521133347722145;    // OSGB false north
            }
            if (grid == "Irish")
            {
                a = 6377340.189;       // OSGBI semi-major
                b = 6356034.447;        // OSGBI semi-minor
                e0 = 200000;           // easting of false origin
                n0 = 250000;          // northing of false origin
                f0 = 1.000035;     // OSGBI scale factor on central meridian
                e2 = 0.00667054015;  // OSGBI eccentricity squared
                lam0 = -0.13962634015954636615389526147909;  // OSGBI false east
                phi0 = 0.93375114981696632365417456114141;    // OSGBI false north
            }

            if (grid == "Channel Islands")
            {
                a = 6378388.000;       // INT24 ED50 semi-major
                b = 6356911.946;       // INT24 ED50 semi-minor 
                e0 = 500000;           // easting of false origin
                n0 = 0;                // northing of false origin
                f0 = 0.9996;           // INT24 ED50 scale factor on central meridian
                e2 = 0.0067226700223333;  // INT24 ED50 eccentricity squared
                lam0 = -0.0523598775598;  // INT24 ED50 false east
                phi0 = 0 * deg2rad;    // INT24 ED50 false north 
            }
            double af0 = a * f0;
            double bf0 = b * f0;
            // easting
            double slat2 = sin(phi) * sin(phi);
            double nu = af0 / (sqrt(1 - (e2 * (slat2))));
            double rho = (nu * (1 - e2)) / (1 - (e2 * slat2));
            double eta2 = (nu / rho) - 1;
            double p = lam - lam0; //LAM-LAM0
            double IV = nu * cos(phi);
            double clat3 = cos(phi) * cos(phi) * cos(phi);
            double tlat2 = tan(phi) * tan(phi);
            double V = (nu / 6) * clat3 * ((nu / rho) - tlat2);
            double clat5 = pow(cos(phi), 5);
            double tlat4 = pow(tan(phi), 4);
            double VI = (nu / 120) * clat5 * ((5 - (18 * tlat2)) + tlat4 + (14 * eta2) - (58 * tlat2 * eta2));
            double east = e0 + (p * IV) + (pow(p, 3) * V) + (pow(p, 5) * VI);
            // northing
            double n = (af0 - bf0) / (af0 + bf0);
            double M = Marc(bf0, n, phi0, phi);
            double I = M + (n0);
            double II = (nu / 2) * sin(phi) * cos(phi);
            double III = ((nu / 24) * sin(phi) * clat3) * (5 - tlat2 + (9 * eta2));
            double IIIA = ((nu / 720) * sin(phi) * clat5) * (61 - (58 * tlat2) + tlat4);
            double north = I + ((p * p) * II) + (pow(p, 4) * III) + (pow(p, 6) * IIIA);
            north = round(north * 100) / 100;
            east = round(east * 100) / 100;

            east = round(east);
            north = round(north);
            EastingValue = east.ToString();
            NorthingValue = north.ToString();
            //en2ngr(east, north, grid); 
        }
        public double Marc(double bf0, double n, double phi0, double phi)
        {
            double Marc = bf0 * (((1 + n + ((5 / 4) * (n * n)) + ((5 / 4) * (n * n * n))) * (phi - phi0))
            - (((3 * n) + (3 * (n * n)) + ((21 / 8) * (n * n * n))) * (sin(phi - phi0)) * (cos(phi + phi0)))
            + ((((15 / 8) * (n * n)) + ((15 / 8) * (n * n * n))) * (sin(2 * (phi - phi0))) * (cos(2 * (phi + phi0))))
            - (((35 / 24) * (n * n * n)) * (sin(3 * (phi - phi0))) * (cos(3 * (phi + phi0)))));
            return (Marc);
        }
        public double abs(double x)
        {
            return Math.Abs(x);
        }
        public double floor(double x)
        {
            return Math.Floor(x);
        }
        public double mod(double y, double x)
        {
            if (y >= 0)
                return y - x * floor(y / x);
            else
                return y + x * (floor(-y / x) + 1.0);
        }
        public double atan2(double y, double x)
        {
            return Math.Atan2(y, x);
        }
        public double sqrt(double x)
        {
            return Math.Sqrt(x);
        }
        public double tan(double x)
        {
            return Math.Tan(x);
        }
        public double sin(double x)
        {
            return Math.Sin(x);
        }
        public double cos(double x)
        {
            return Math.Cos(x);
        }
        public double acos(double x)
        {
            return Math.Acos(x);
        }
        public double round(double x)
        {
            return Math.Round(x);
        }
        public double ceil(double x)
        {
            return Math.Ceiling(x);
        }
        public double ln(double x)
        {
            return Math.Log(x);
        }
        public double pow(double x, double y)
        {
            return Math.Pow(x, y);
        }
        public double atan(double x)
        {
            return Math.Atan(x);
        }
        public double Number(String x)
        {
            return double.Parse(x);
        }

    }
}

/**
 * Geoid.java
 * A Java implementation of the GeographicLib::Geoid class.
 *
 * <p>
 * This file is a translation of GeographicLib/src/Geoid.cpp from C++ to Java.
 * Original C++ implementation by Charles Karney (2009-2020) <karney@alum.mit.edu>
 * and licensed under the MIT/X11 License. For more information, see:
 * <a href="https://geographiclib.sourceforge.io/">GeographicLib Project</a>.
 * </p>
 *
 * <p>
 * Ported to Java by [Akın BÜYÜKBULUT] (<akinbuyukbulut@gmail.com>).
 * This work is licensed under the MIT License.
 * </p>
 *
 * <p>
 * This file implements methods for geoid height calculations, caching,
 * and interpolation, providing functionality equivalent to the original C++
 * Geoid class in the GeographicLib library.
 * </p>
 *
 * @author [Akın BÜYÜKBULUT]
 * @version 1.0
 * @since 2025
 */

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static java.lang.Math.*;

public final class Geoid {
    private static final int pixel_size_ = 2;
    private static final long pixel_max_ = 0xFFFFL;  // 65.535 as unsigned decimal

    /*
     #if GEOGRAPHICLIB_GEOID_PGM_PIXEL_WIDTH != 4
     pixel_size_ = 4;
     pixel_max_ = 0xFFFFFFFFL; // 4.294.967.295 as unsigned decimal
     */

    private static final int stencilsize_ = 12;
    private static final int nterms_ = ((3 + 1) * (3 + 2)) / 2; // for a cubic fit

    private static final int c0_ = 240;
    private static final int c0n_ = 372;
    private static final int c0s_ = 372;

    private static final double TD = 360.0; //Total Degree
    private static final double QD = 90.0; //Quarter Degree
    private static final double HD = 180.0; //Half Degree
    private static final int[] c3_ = {
            9, -18, -88,    0,  96,   90,   0,   0, -60, -20,
            -9,  18,   8,    0, -96,   30,   0,   0,  60, -20,
            9, -88, -18,   90,  96,    0, -20, -60,   0,   0,
            186, -42, -42, -150, -96, -150,  60,  60,  60,  60,
            54, 162, -78,   30, -24,  -90, -60,  60, -60,  60,
            -9, -32,  18,   30,  24,    0,  20, -60,   0,   0,
            -9,   8,  18,   30, -96,    0, -20,  60,   0,   0,
            54, -78, 162,  -90, -24,   30,  60, -60,  60, -60,
            -54,  78,  78,   90, 144,   90, -60, -60, -60, -60,
            9,  -8, -18,  -30, -24,    0,  20,  60,   0,   0,
            -9,  18, -32,    0,  24,   30,   0,   0, -60,  20,
            9, -18,  -8,    0, -24,  -30,   0,   0,  60,  20,
    };

    private static final int[] c3n_ = {
            0, 0, -131, 0,  138,  144, 0,   0, -102, -31,
            0, 0,    7, 0, -138,   42, 0,   0,  102, -31,
            62, 0,  -31, 0,    0,  -62, 0,   0,    0,  31,
            124, 0,  -62, 0,    0, -124, 0,   0,    0,  62,
            124, 0,  -62, 0,    0, -124, 0,   0,    0,  62,
            62, 0,  -31, 0,    0,  -62, 0,   0,    0,  31,
            0, 0,   45, 0, -183,   -9, 0,  93,   18,   0,
            0, 0,  216, 0,   33,   87, 0, -93,   12, -93,
            0, 0,  156, 0,  153,   99, 0, -93,  -12, -93,
            0, 0,  -45, 0,   -3,    9, 0,  93,  -18,   0,
            0, 0,  -55, 0,   48,   42, 0,   0,  -84,  31,
            0, 0,   -7, 0,  -48,  -42, 0,   0,   84,  31,
    };

    private static final int[] c3s_ = {
            18,  -36, -122,   0,  120,  135,   0,   0,  -84, -31,
            -18,   36,   -2,   0, -120,   51,   0,   0,   84, -31,
            36, -165,  -27,  93,  147,   -9,   0, -93,   18,   0,
            210,   45, -111, -93,  -57, -192,   0,  93,   12,  93,
            162,  141,  -75, -93, -129, -180,   0,  93,  -12,  93,
            -36,  -21,   27,  93,   39,    9,   0, -93,  -18,   0,
            0,    0,   62,   0,    0,   31,   0,   0,    0, -31,
            0,    0,  124,   0,    0,   62,   0,   0,    0, -62,
            0,    0,  124,   0,    0,   62,   0,   0,    0, -62,
            0,    0,   62,   0,    0,   31,   0,   0,    0, -31,
            -18,   36,  -64,   0,   66,   51,   0,   0, -102,  31,
            18,  -36,    2,   0,  -66,  -51,   0,   0,  102,  31,
    };

    // Flags indicating conversions between heights above the geoid and heights above the ellipsoid.
    enum ConvertFlag {
        ELLIPSOIDTOGEOID(-1), //The multiplier for converting from heights above the ellipsoid to heights above the geoid.
        NONE(0), //No conversion.
        GEOIDTOELLIPSOID(1); //The multiplier for converting from heights above the geoid to heights above the ellipsoid.

        private final int value;

        ConvertFlag(int v) {
            value = v;
        }

        public int getValue() {
            return value;
        }
    }

    private final String _name;
    private String _dir;
    private final String _filename;
    private final double EQUATORIAL_RADIUS = 6378137.0;
    private final double FLATTENING = 1.0 / 298.257223563;
    private final boolean _cubic;
    private RandomAccessFile _file;
    private final double _rlonres, _rlatres;
    private final String _description;
    private final String _datetime;
    private final double _offset, _scale, _maxerror, _rmserror;
    private final int _width, _height;
    private final long _datastart, _swidth;
    private boolean _threadsafe;
    private final List<int[]> _data = new ArrayList<>();
    private boolean _cache;
    private int _xoffset, _yoffset, _xsize, _ysize;
    private int _ix, _iy;
    private double _v00, _v01, _v10, _v11;
    private final double[] _t = new double[nterms_];
    public static final String DEFAULT_GEOID_PATH;
    public static final String DEFAULT_GEOID_NAME;

    static {
        final String GEOGRAPHICLIB_GEOID_DEFAULT_NAME = "egm2008-5";
        DEFAULT_GEOID_PATH = getGeoidPath();
        DEFAULT_GEOID_NAME = getGeoidName(GEOGRAPHICLIB_GEOID_DEFAULT_NAME);
    }

    /**
     * Constructs a Geoid object.
     *
     * @param name       The name of the geoid.
     * @param path       (Optional) The directory for the data file. If not specified or empty,
     *                   the default path will be used.
     * @param cubic      The interpolation method. If {@code false}, bilinear interpolation
     *                   is used. If {@code true}, cubic interpolation is used.
     * @param threadsafe If {@code true}, creates a thread-safe object.
     * @throws RuntimeException If the data file cannot be found, is unreadable, or is corrupt.
     *
     *                          <p>The data file is identified by appending ".pgm" to the {@code name}. If {@code path} is
     *                          specified and non-empty, the file is loaded from the specified directory. Otherwise, the
     *                          default path is used as determined by {@code DefaultGeoidPath()}.</p>
     *
     *                          <p>If {@code threadsafe} is {@code true}, the dataset is fully loaded into memory, the data
     *                          file is closed, and single-cell caching is turned off, resulting in a thread-safe {@code Geoid}
     *                          object.</p>
     */
    public Geoid(String name, String path, boolean cubic, boolean threadsafe) throws IOException {
        _cubic = cubic;
        _threadsafe = threadsafe;
        _name = name;
        _dir = path;
        if (_dir == null || _dir.trim().isEmpty()) {
            _dir = DEFAULT_GEOID_PATH;
        }
        _filename = _dir + File.separator + _name + ".pgm";
        _file = new RandomAccessFile(_filename, "r");
        int dataWidth = 0, dataHeight = 0;
        double offset = Double.MAX_VALUE, scale = 0, maxError = -1, rmsError = -1;
        String description = "NONE";
        String dateTime = "UNKNOWN";
        _file.seek(0);
        byte[] headerBytes = new byte[2048];
        int bytesRead = _file.read(headerBytes);
        if (bytesRead <= 0) {
            throw new RuntimeException("File is empty or cannot be read: " + _filename);
        }

        String headerText = new String(headerBytes, 0, bytesRead, StandardCharsets.UTF_8);
        BufferedReader bufferedReader = new BufferedReader(new java.io.StringReader(headerText));
        String line = bufferedReader.readLine();
        if (!"P5".equals(line)) {
            throw new RuntimeException("File not in PGM format: " + _filename);
        }

        boolean rasterSizeRead = false;
        boolean maxValueRead = false;

        // in ASCII decimal. Must be less than 65536 and more than zero.
        long maximumGrayValue = 0L;

        Pattern pattern = Pattern.compile("^#\\s+([A-Za-z]+)\\s+(.+)$");
        while ((line = bufferedReader.readLine()) != null) {
            if (line.startsWith("#")) {
                Matcher matcher = pattern.matcher(line);
                if (matcher.find()) {
                    String key = matcher.group(1);
                    String value = matcher.group(2).trim();
                    switch (key) {
                        case "Description":
                            description = value;
                            break;
                        case "DateTime":
                            dateTime = value;
                            break;
                        case "Offset": {
                            try {
                                offset = Double.parseDouble(value);
                            } catch (Exception e) {
                                throw new RuntimeException("Error reading offset: " + _filename);
                            }
                            break;
                        }
                        case "Scale": {
                            try {
                                scale = Double.parseDouble(value);
                            } catch (Exception e) {
                                throw new RuntimeException("Error reading scale: " + _filename);
                            }
                            break;
                        }
                        case "MaxCubicError":
                        case "MaxBilinearError": {
                            if ((_cubic && key.equals("MaxCubicError")) || (!_cubic && key.equals("MaxBilinearError"))) {
                                try {
                                    maxError = Double.parseDouble(value);
                                } catch (Exception e) {
                                    System.out.println("Warning: Error reading " + key + ": " + _filename);
                                }
                            }
                            break;
                        }
                        case "RMSCubicError":
                        case "RMSBilinearError": {
                            if ((_cubic && key.equals("RMSCubicError")) || (!_cubic && key.equals("RMSBilinearError"))) {
                                try {
                                    rmsError = Double.parseDouble(value);
                                } catch (Exception e) {
                                    System.out.println("Warning: Error reading " + key + ": " + _filename);
                                }
                            }
                            break;
                        }
                        default:
                            // TODO: Handle other cases if needed
                            // URL, ARE_OR_POINT, Vertical_Datum etc.
                            break;
                    }
                }
            } else {
                if (!rasterSizeRead) {
                    String[] items = line.split("\\s+");
                    if (items.length != 2) {
                        throw new RuntimeException("Error reading raster size: " + _filename);
                    }
                    try {
                        dataWidth = Integer.parseInt(items[0]);
                    } catch (NumberFormatException e) {
                        throw new RuntimeException("Error reading raster width: " + _filename);
                    }
                    try {
                        dataHeight = Integer.parseInt(items[1]);
                    } catch (NumberFormatException e) {
                        throw new RuntimeException("Error reading raster height: " + _filename);
                    }
                    rasterSizeRead = true;
                } else {
                    try {
                        maximumGrayValue = Long.parseLong(line);
                    } catch (NumberFormatException e) {
                        throw new RuntimeException("Error reading max value: " + _filename);
                    }
                    maxValueRead = true;
                    break;
                }
            }
        }

        if (!maxValueRead) {
            throw new RuntimeException("Could not read maximum gray value in PGM header: " + _filename);
        }

        if (maximumGrayValue != pixel_max_) {
            throw new RuntimeException("Incorrect value of maximum gray value " + _filename);
        }

        // Add 1 for whitespace after maxval
        long dataStart = headerText.indexOf(String.valueOf(maximumGrayValue)) + String.valueOf(maximumGrayValue).length() + 1;
        _file.seek(dataStart);

        if (offset == Double.MAX_VALUE) {
            throw new RuntimeException("Offset not set: " + _filename);
        }
        if (scale == 0) {
            throw new RuntimeException("Scale not set: " + _filename);
        }
        if (scale < 0) {
            throw new RuntimeException("Scale must be positive: " + _filename);
        }
        if (dataHeight < 2 || dataWidth < 2) {
            throw new RuntimeException("Raster size too small: " + _filename);
        }
        if ((dataWidth & 1) != 0) {
            throw new RuntimeException("Raster width is odd: " + _filename);
        }
        if ((dataHeight & 1) == 0) {
            throw new RuntimeException("Raster height is even: " + _filename);
        }


        long fileLength = _file.length();
        long expectedDataLength = dataStart + pixel_size_ * ((long) dataWidth) * ((long) dataHeight);
        if (expectedDataLength != fileLength) {
            throw new RuntimeException("File has the wrong length: " + _filename);
        }

        _offset = offset;
        _scale = scale;
        _maxerror = maxError;
        _rmserror = rmsError;
        _description = description;
        _datetime = dateTime;
        _width = dataWidth;
        _height = dataHeight;
        _datastart = dataStart;
        _swidth = dataWidth;
        _rlonres = _width / TD;
        _rlatres = (_height - 1) / HD;
        _cache = false;
        _ix = _width;
        _iy = _height;

        if (threadsafe) {
            cacheAll();
            _file.close();
            _threadsafe = true;
        } else {
            _threadsafe = false;
        }
    }

    /**
     * Sets up a cache for the specified rectangular area.
     *
     * @param south The latitude (in degrees) of the southern edge of the cached area.
     *              Must be in the range [-90°, 90°].
     * @param west  The longitude (in degrees) of the western edge of the cached area.
     * @param north The latitude (in degrees) of the northern edge of the cached area.
     *              Must be in the range [-90°, 90°].
     * @param east  The longitude (in degrees) of the eastern edge of the cached area.
     *              Interpreted as being east of {@code west}, adding 360° to its value if necessary.
     * @throws RuntimeException If there is a problem reading the data.
     * @throws RuntimeException If this method is called on a thread-safe Geoid.
     *
     *                          <p>Caches the data for the specified "rectangular" area bounded by the parallels
     *                          {@code south} and {@code north}, and the meridians {@code west} and {@code east}.
     *                          </p>
     *
     *                          <p>The {@code east} parameter is always interpreted as being east of {@code west}. If necessary,
     *                          360° is added to {@code east} to ensure this condition. The latitude parameters {@code south}
     *                          and {@code north} must lie within the range [-90°, 90°].</p>
     */
    public void cacheArea(double south, double west, double north, double east) {
        if (_threadsafe) {
            throw new RuntimeException("Attempt to change cache of threadsafe Geoid");
        }
        if (south > north) {
            cacheClear();
            return;
        }
        south = clampLatitude(south);
        north = clampLatitude(north);
        west = normalizeLongitude(west);
        east = normalizeLongitude(east);
        if (east <= west) {
            east += TD;
        }
        int indexWest = (int) floor(west * _rlonres);
        int indexEast = (int) floor(east * _rlonres);
        int indexNorth = (int) floor(-north * _rlatres) + (_height - 1) / 2;
        int indexSouth = (int) floor(-south * _rlatres) + (_height - 1) / 2;

        indexNorth = max(0, min(_height - 2, indexNorth));
        indexSouth = max(0, min(_height - 2, indexSouth));
        indexSouth += 1;
        indexEast += 1;

        if (_cubic) {
            indexNorth -= 1;
            indexSouth += 1;
            indexWest -= 1;
            indexEast += 1;
        }
        if (indexEast - indexWest >= _width - 1) {
            indexWest = 0;
            indexEast = _width - 1;
        } else {
            if (indexWest < 0) {
                indexWest += _width;
            } else if (indexWest >= _width) {
                indexWest -= _width;
            }
            if (indexEast < 0) {
                indexEast += _width;
            } else if (indexEast >= _width) {
                indexEast -= _width;
            }
        }

        _xsize = indexEast - indexWest + 1;
        _ysize = indexSouth - indexNorth + 1;
        _xoffset = indexWest;
        _yoffset = indexNorth;

        _data.clear();
        for (int iy = 0; iy < _ysize; iy++) {
            _data.add(new int[_xsize]);
        }

        try {
            for (int iy = indexNorth; iy <= indexSouth; ++iy) {
                int iy1 = iy;
                int iw1 = indexWest;
                if (iy1 < 0 || iy1 >= _height) {
                    iy1 = (iy1 < 0) ? -iy1 : 2 * (_height - 1) - iy1;
                    iw1 += _width / 2;
                    if (iw1 >= _width) {
                        iw1 -= _width;
                    }
                }
                int xs1 = min(_width - iw1, _xsize);
                filepos(iw1, iy1);
                int[] rowData = _data.get(iy - indexNorth);
                readShortsToIntArray(_file, rowData, 0, xs1);
                if (xs1 < _xsize) {
                    filepos(0, iy1);
                    readShortsToIntArray(_file, rowData, xs1, _xsize - xs1);
                }
            }
            _cache = true;
        } catch (IOException e) {
            cacheClear();
            throw new RuntimeException("Error filling cache: " + e.getMessage());
        }
    }

    /**
     * Caches all the data for the geoid model.
     *
     * @throws RuntimeException If there is a problem reading the data.
     * @throws RuntimeException If this method is called on a thread-safe Geoid.
     *
     *                          <p>This method is efficient for data sets with a grid resolution of 5' or coarser. Note that
     *                          since the cache values are stored as integers in Java (instead of unsigned short in the original
     *                          implementation), the memory requirements are approximately double compared to the original
     *                          GeographicLib implementation:</p>
     *                          <ul>
     *                            <li>1' grid: ~900MB</li>
     *                            <li>2.5' grid: ~144MB</li>
     *                            <li>5' grid: ~36MB</li>
     *                          </ul>
     *                          <p>Consider this when working with large grids to ensure sufficient memory is available.</p>
     */
    public void cacheAll() {
        cacheArea(-QD, 0, QD, TD);
    }

    /**
     * Clears the cache for the geoid model.
     *
     * <p>This method never throws an error. If called on a thread-safe Geoid, this method does nothing.</p>
     */
    public void cacheClear() {
        if (!_threadsafe) {
            _cache = false;
            _data.clear();
        }
    }

    /**
     * Converts a height above the geoid to a height above the ellipsoid, or vice versa.
     *
     * @param latitude The latitude of the point (degrees).
     * @param longitude The longitude of the point (degrees).
     * @param heightOfThePoint The height of the point (meters).
     * @param direction A {@code ConvertFlag} specifying the direction of the conversion:
     *                  {@code GEOID_TO_ELLIPSOID} converts a height above the geoid to a height
     *                  above the ellipsoid, while {@code ELLIPSOID_TO_GEOID} converts a height
     *                  above the ellipsoid to a height above the geoid.
     * @return The converted height (meters).
     */
    public double convertHeight(double latitude, double longitude, double heightOfThePoint, ConvertFlag direction) {
        return heightOfThePoint + direction.getValue() * height(latitude, longitude);
    }

    /**
     * Computes the height of the geoid above the ellipsoid at a specified point.
     *
     * @param latitude The latitude of the point (degrees). Must be in the range [-90°, 90°].
     * @param longitude The longitude of the point (degrees).
     * @return The height of the geoid above the ellipsoid (meters).
     */
    public double computeGeoidHeight(double latitude, double longitude) {
        return height(latitude, longitude);
    }

    public String getDescription() {
        return _description;
    }

    public String getDateTime() {
        return _datetime;
    }

    public String getGeoidName() {
        return _name;
    }

    public String getGeoidFile() {
        return _filename;
    }

    public String getGeoidDirectory() {
        return _dir;
    }

    public String getInterpolation() {
        return _cubic ? "cubic" : "bilinear";
    }

    public double getMaxError() {
        return _maxerror;
    }

    public double getRMSError() {
        return _rmserror;
    }

    public double getOffset() {
        return _offset;
    }

    public double getScale() {
        return _scale;
    }

    public boolean isThreadSafe() {
        return _threadsafe;
    }

    public boolean isCacheEnabled() {
        return _cache;
    }

    public double getEquatorialRadius() {
        return EQUATORIAL_RADIUS;
    }

    public double getFlattening() {
        return FLATTENING;
    }

    public double getCacheWest() {
        if (!_cache) {
            return 0;
        }
        return (((_xoffset + (_xsize == _width ? 0 : (_cubic ? 1 : 0)) + _width / 2) % _width) - _width / 2) / _rlonres;
    }

    public double getCacheEast() {
        if (!_cache) {
            return 0;
        }
        return getCacheWest() + (_xsize - (_xsize == _width ? 0 : 1 + 2 * (_cubic ? 1 : 0))) / _rlonres;
    }

    public double getCacheNorth() {
        if (!_cache) {
            return 0;
        }
        return QD - (_yoffset + (_cubic ? 1 : 0)) / _rlatres;
    }

    public double getCacheSouth() {
        if (!_cache) {
            return 0;
        }
        return QD - (_yoffset + _ysize - 1 - (_cubic ? 1 : 0)) / _rlatres;
    }


    /* PRIVATE METHODS */

    private void filepos(int ix, int iy) {
        try {
            long pos = _datastart + pixel_size_ * ((long) iy * _swidth + (long) ix);
            _file.seek(pos);
        } catch (IOException e) {
            throw new RuntimeException("Error seeking file: " + e.getMessage());
        }
    }

    private double rawval(int ix, int iy) {
        if (ix < 0) ix += _width;
        else if (ix >= _width) ix -= _width;

        if (_cache && iy >= _yoffset && iy < _yoffset + _ysize) {
            int localIx = ix;
            if (localIx < _xoffset) {
                localIx += _width;
            }
            localIx -= _xoffset;
            if (localIx >= 0 && localIx < _xsize) {
                return _data.get(iy - _yoffset)[localIx];
            }
        }

        if (iy < 0 || iy >= _height) {
            iy = (iy < 0) ? -iy : 2 * (_height - 1) - iy;
            ix += (ix < _width / 2 ? 1 : -1) * (_width / 2);
        }
        try {
            filepos(ix, iy);
            int hi = _file.read();
            int lo = _file.read();
            if (hi < 0 || lo < 0) {
                throw new RuntimeException("Error reading " + _filename + ": unexpected EOF");
            }
            return ((hi & 0xFF) << 8) | (lo & 0xFF);
        } catch (IOException e) {
            throw new RuntimeException("Error reading " + _filename + ": " + e.getMessage());
        }
    }

    private double height(double latitude, double longitude) {
        latitude = clampLatitude(latitude);
        if (Double.isNaN(latitude) || Double.isNaN(longitude)) {
            return Double.NaN;
        }
        longitude = normalizeLongitude(longitude);

        double fractionalX = longitude * _rlonres;
        double fractionalY = -latitude * _rlatres;
        int ix = (int) floor(fractionalX);
        int iy = min((_height - 1) / 2 - 1, (int) floor(fractionalY));
        fractionalX -= ix;
        fractionalY -= iy;
        iy += (_height - 1) / 2;

        if (ix < 0) ix += _width;
        else if (ix >= _width) ix -= _width;

        double v00 = 0, v01 = 0, v10 = 0, v11 = 0;
        double[] localT = new double[nterms_];

        boolean sameCell = !_threadsafe && (ix == _ix && iy == _iy);

        if (!_cubic) {
            if (!sameCell) {
                v00 = rawval(ix, iy);
                v01 = rawval(ix + 1, iy);
                v10 = rawval(ix, iy + 1);
                v11 = rawval(ix + 1, iy + 1);
            } else {
                v00 = _v00;
                v01 = _v01;
                v10 = _v10;
                v11 = _v11;
            }
        } else {
            if (!sameCell) {
                double[] v = new double[stencilsize_];
                int k = 0;
                v[k++] = rawval(ix, iy - 1);
                v[k++] = rawval(ix + 1, iy - 1);
                v[k++] = rawval(ix - 1, iy);
                v[k++] = rawval(ix, iy);
                v[k++] = rawval(ix + 1, iy);
                v[k++] = rawval(ix + 2, iy);
                v[k++] = rawval(ix - 1, iy + 1);
                v[k++] = rawval(ix, iy + 1);
                v[k++] = rawval(ix + 1, iy + 1);
                v[k++] = rawval(ix + 2, iy + 1);
                v[k++] = rawval(ix, iy + 2);
                v[k++] = rawval(ix + 1, iy + 2);

                int[] c3x;
                int c0x;
                if (iy == 0) {
                    c3x = c3n_;
                    c0x = c0n_;
                } else if (iy == _height - 2) {
                    c3x = c3s_;
                    c0x = c0s_;
                } else {
                    c3x = c3_;
                    c0x = c0_;
                }

                for (int i = 0; i < nterms_; i++) {
                    double sum = 0;
                    for (int j = 0; j < stencilsize_; j++) {
                        sum += v[j] * c3x[nterms_ * j + i];
                    }
                    localT[i] = sum / c0x;
                }
            } else {
                System.arraycopy(_t, 0, localT, 0, nterms_);
            }
        }

        double hVal;
        if (!_cubic) {
            double a = (1 - fractionalX) * v00 + fractionalX * v01;
            double b = (1 - fractionalX) * v10 + fractionalX * v11;
            double c = (1 - fractionalY) * a + fractionalY * b;
            hVal = _offset + _scale * c;
            if (!sameCell) {
                _ix = ix;
                _iy = iy;
                _v00 = v00;
                _v01 = v01;
                _v10 = v10;
                _v11 = v11;
            }
        } else {
            double h = localT[0] + fractionalX * (localT[1] + fractionalX * (localT[3] + fractionalX * localT[6])) + fractionalY * (localT[2] + fractionalX * (localT[4] + fractionalX * localT[7]) + fractionalY * (localT[5] + fractionalX * localT[8] + fractionalY * localT[9]));
            hVal = _offset + _scale * h;
            if (!sameCell) {
                _ix = ix;
                _iy = iy;
                System.arraycopy(localT, 0, _t, 0, nterms_);
            }
        }
        return hVal;
    }

    private static void readShortsToIntArray(RandomAccessFile randomAccessFile, int[] destination, int offset, int length) throws IOException {
        for (int index = 0; index < length; index++) {
            int high = randomAccessFile.read();
            int low = randomAccessFile.read();
            if (high < 0 || low < 0) {
                throw new RuntimeException("Error reading file - unexpected EOF");
            }
            int val = ((high & 0xFF) << 8) | (low & 0xFF);
            destination[offset + index] = val;
        }
    }


    private static double normalizeLongitude(double longitude) {

        double y = Math.IEEEremainder(longitude, TD);
        return Math.abs(y) == HD ? Math.copySign(HD, longitude) : y;
    }

    private static double clampLatitude(double latitude) {
        return max(-90.0, min(90.0, latitude));
    }

    private static String getGeoidPath() {
        String path = System.getenv("GEOGRAPHICLIB_GEOID_PATH");
        if (path != null && !path.trim().isEmpty()) {
            return path;
        }
        path = System.getenv("GEOGRAPHICLIB_DATA");
        if (path != null && !path.trim().isEmpty()) {
            return path + File.separator + "geoids";
        }
        return "geoids";
    }

    private static String getGeoidName(String defaultName) {
        String name = System.getenv("GEOGRAPHICLIB_GEOID_NAME");
        return (name != null && !name.trim().isEmpty()) ? name : defaultName;
    }
}

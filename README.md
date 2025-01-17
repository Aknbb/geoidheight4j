GeoidHeight4J
==============

**GeoidHeight4J** is a pure Java library that provides tools to compute geoid heights
and convert heights between the geoid and ellipsoid. This library is a Java port
of the Geoid class from the [GeographicLib project.](https://geographiclib.sourceforge.io/) **GeoidHeight4J** is ideal for applications requiring precise geoid height calculations,
such as geographic information systems (GIS) and geodetic applications.

The reason why the class variables in this class are written with abbreviations like _rlonres and _v00 is that they are intentionally kept the same as the variable names in the [Geoid.hpp](https://github.com/ObjSal/GeographicLib/blob/master/include/GeographicLib/Geoid.hpp) class.

You can use the [Geoid Height Document](https://geographiclib.sourceforge.io/C++/doc/geoid.html) and [Geoid Class Reference](https://geographiclib.sourceforge.io/C++/doc/classGeographicLib_1_1Geoid.html) links to access EGM data and learn more about the Geoid.
 
Features
--------

- Geoid Height Calculation: Compute the height of the geoid above the ellipsoid
  for a given latitude and longitude.
- Height Conversion: Convert heights between the geoid and ellipsoid reference
  surfaces.
- Interpolation Methods: Supports bilinear and cubic interpolation for height
  computations.

Usage
-----

Example Code:

```java
import your.package.Geoid;

public class Main {
    public static void main(String[] args) {
        // Initialize the Geoid object with the model name and path to data
        Geoid geoid = new Geoid("egm2008-5", "/home/aknbb/Desktop/geoid-data", true, false);

        // Example: Compute the geoid height at a specific latitude and longitude
        double latitude = 37.97150021753285;  // Latitude in degrees
        double longitude = 34.62309350397177;  // Longitude in degrees
        double geoidHeight = geoid.computeGeoidHeight(latitude, longitude);
        // 35.20602975851298 meters for egm2008-5.pgm data with cubic interpolation
        System.out.println("Geoid height above ellipsoid: " + geoidHeight + " meters"); 
    }
}
```

API Overview
------------

| Method                                                                        | Description                                                                                |
|-------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------|
| Geoid(String name, String path, boolean cubic, boolean threadsafe)            | Initializes the Geoid object with the given parameters.                                    |
| double computeGeoidHeight(double latitude, double longitude)                             | Computes the height of the geoid above the ellipsoid for the given latitude and longitude. |
| double convertHeight(double latitude, double longitude, double heightOfThePoint, ConvertFlag direction) | Converts height between geoid and ellipsoid references.                                    |
| void cacheArea(double south, double west, double north, double east)          | Caches the geoid data for a specific area to improve performance.                          |
| void cacheAll()                                                               | Caches all geoid data into memory to optimize computations across the entire dataset.    |
| void cacheClear()                                                             | Clears the cache for the geoid model.    |

For more detailed information, refer to the Javadoc content in the Geoid class or visit [GeographicLib's documentation.](https://geographiclib.sourceforge.io/C++/doc/classGeographicLib_1_1Geoid.html)

### License

This project is licensed under the [MIT License](LICENSE).

#### Contact
Feel free to explore my Tile Layout Converter project and get in touch if you have any questions or collaboration ideas. You can reach out to me via [akinbuyukbulut@gmail.com](mailto:akinbuyukbulut@gmail.com) or connect with me on [LinkedIn](https://www.linkedin.com/in/akinbuyukbulut/) and [GitHub](https://github.com/Aknbb).
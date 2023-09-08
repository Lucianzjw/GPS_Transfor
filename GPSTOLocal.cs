using System;
using System.Collections;
using UnityEngine;

public struct Vector2D
{
    public double x;
    public double y;

    public Vector2D(double xValue, double yValue)
    {
        x = xValue;
        y = yValue;
    }
}

public class GPSTOLocal 
{

    #region 坐标

    public Vector2D BottomRightPoint = new Vector2D(240.2, -156.2); //Unity中右下点  （X正方向和Z轴的负方向之间）
    //public Vector2D TopLeftPoint = new Vector2D(380, -85.4); //Unity中左上点  （Z轴正方向和X轴负方向之间）
    public Vector2D TopLeftPoint = new Vector2D(-239.4,82.7); //Unity中左上点  （Z轴正方向和X轴负方向之间）

    public Vector2D BottomRightSai = new Vector2D(117.127171262493, 31.8290517125821); //地图中对应的左上经纬度点  
    //public Vector2D TopLeftSai = new Vector2D(117.12572873735844, 31.82840269630183); //地图中对应的右下经纬度点  
    public Vector2D TopLeftSai = new Vector2D( 117.13223928613924, 31.826901967189368); //地图中对应的右下经纬度点  

    private double z_offset, x_offset, z_w_offset, x_w_offset;

    protected override void Awake()
    {
        base.Awake();
        InitBasicNum();
    }
    

    private void InitBasicNum()
    {

        z_offset = BottomRightSai.y - TopLeftSai.y; //地图中的维度差  
        x_offset = BottomRightSai.x - TopLeftSai.x; //地图中的经度差  
        z_w_offset = BottomRightPoint.y - TopLeftPoint.y; //unity中的维度差  
        x_w_offset = BottomRightPoint.x - TopLeftPoint.x; //unity中的经度差  
    }

    /// <summary>
    /// 由经纬度得到位置点  
    /// </summary>
    /// <param name="se"></param>
    /// <returns></returns>
    public Vector2D GetWorldPoint(Vector2D se)
    {
        double tempX = se.y - TopLeftSai.x;
        double tempZ = se.x - BottomRightSai.y;
        double _tempX = (tempX * x_w_offset / x_offset + TopLeftPoint.x);
        double _tempZ = (tempZ * z_w_offset / z_offset + BottomRightPoint.y);
        //坐标偏差（在Unity中的坐标）
        return new Vector2D(_tempX, _tempZ);
    }

    /// <summary>
    /// 由位置点得到经纬度  
    /// </summary>
    /// <param name="curPoint"></param>
    /// <returns></returns>
    public Vector2D GetLatLon(Vector2D curPoint)
    {
        //坐标偏差
        double _x_offset = (curPoint.x - BottomRightPoint.x) * x_offset / x_w_offset;
        double _z_offset = (curPoint.y - TopLeftPoint.y) * z_offset / z_w_offset;
        double resultX = _x_offset + BottomRightSai.x;
        double resultZ = _z_offset + TopLeftSai.y;
        return new Vector2D(resultZ,resultX);
    }

    #endregion

    #region 墨卡托


    /// <summary>
    /// 墨卡托坐标转换（经纬度与米之间互转）
    /// </summary>

    /// <summary>
    /// 经纬度转Web墨卡托（单位：米）
    /// </summary>
    /// <param name="longitude">经度</param>
    /// <param name="latitude">纬度</param>
    /// <returns>转换后的位置</returns>
    public Vector2D Degree2WebMercatorMeter(double latitude, double longitude)
    {
        var xValue = longitude * 20037508.34 / 180;
        var y = Math.Log(Math.Tan((90 + latitude) * Math.PI / 360)) / (Math.PI / 180);
        var yValue = y * 20037508.34 / 180;
        xValue -= 13030000;
        yValue -= 3741000;
        return new Vector2D(xValue, yValue);
    }


    /// <summary>
    /// 经纬度转World墨卡托（单位：米）
    /// </summary>
    /// <param name="longitude">经度</param>
    /// <param name="latitude">纬度</param>
    /// <returns>转换后的位置</returns>
    public Vector2D Degree2WorldMercatorMeter(double longitude, double latitude)
    {
        const int radius = 6378137;
        const double minorRadius = 6356752.314245179;

        const double d = Math.PI / 180;
        const double r = radius;
        var y = latitude * d;
        const double tmp = minorRadius / r;
        double e = Math.Sqrt(1 - tmp * tmp),
            con = e * Math.Sin(y);

        var ts = Math.Tan(Math.PI / 4 - y / 2) / Math.Pow((1 - con) / (1 + con), e / 2);
        y = -r * Math.Log(Math.Max(ts, 1E-10));

        var xValue = longitude * d * r;
        var yValue = y;

        return new Vector2D(xValue, yValue);
    }

    /// <summary>
    /// Web墨卡托转经纬度
    /// </summary>
    /// <param name="x">X坐标值（单位：米）</param>
    /// <param name="y">Y坐标值（单位：米）</param>
    /// <returns>转换后的位置</returns>
    public Vector2D WebMercatorMeter2Degree(double x, double y)
    {
        var xValue = x / 20037508.34 * 180;
        var yValue = y / 20037508.34 * 180;
        yValue = 180 / Math.PI * (2 * Math.Atan(Math.Exp(yValue * Math.PI / 180)) - Math.PI / 2);
        var longitude = xValue;
        var latitude = yValue;
        return new Vector2D(longitude, latitude);
    }


    #endregion

    #region WGS2GCJ

    private const double pi = 3.1415926535897932384626;
    private const double x_pi = 3.14159265358979324 * 3000.0 / 180.0;
    private const double a = 6378245.0;
    private const double ee = 0.00669342162296594323;

    public double transformLat(double x, double y)
    {
        double ret = -100.0 + 2.0 * x + 3.0 * y + 0.2 * y * y + 0.1 * x * y
                     + 0.2 * Math.Sqrt(Math.Abs(x));
        ret += (20.0 * Math.Sin(6.0 * x * pi) + 20.0 * Math.Sin(2.0 * x * pi)) * 2.0 / 3.0;
        ret += (20.0 * Math.Sin(y * pi) + 40.0 * Math.Sin(y / 3.0 * pi)) * 2.0 / 3.0;
        ret += (160.0 * Math.Sin(y / 12.0 * pi) + 320 * Math.Sin(y * pi / 30.0)) * 2.0 / 3.0;
        return ret;
    }

    public double transformLon(double x, double y)
    {
        double ret = 300.0 + x + 2.0 * y + 0.1 * x * x + 0.1 * x * y + 0.1
            * Math.Sqrt(Math.Abs(x));
        ret += (20.0 * Math.Sin(6.0 * x * pi) + 20.0 * Math.Sin(2.0 * x * pi)) * 2.0 / 3.0;
        ret += (20.0 * Math.Sin(x * pi) + 40.0 * Math.Sin(x / 3.0 * pi)) * 2.0 / 3.0;
        ret += (150.0 * Math.Sin(x / 12.0 * pi) + 300.0 * Math.Sin(x / 30.0
                                                                   * pi)) * 2.0 / 3.0;
        return ret;
    }

    /// <summary>
    /// 判断当前坐标是否是在中国内
    /// </summary>
    /// <param name="lat"></param>
    /// <param name="lon"></param>
    /// <returns></returns>
    public bool outOfChina(double lat, double lon)
    {
        if (lon < 72.004 || lon > 137.8347)
            return true;
        if (lat < 0.8293 || lat > 55.8271)
            return true;
        return false;
    }

    /// <summary>
    /// WGS84 to 火星坐标系 (GCJ-02) 
    /// </summary>
    /// <param name="lat"></param>
    /// <param name="lon"></param>
    /// <returns></returns>
    public Vector2D gps84_To_Gcj02(double lat, double lon)
    {
        if (outOfChina(lat, lon))
        {
            return new Vector2D(lat, lon);
        }

        double dLat = transformLat(lon - 105.0, lat - 35.0);
        double dLon = transformLon(lon - 105.0, lat - 35.0);
        double radLat = lat / 180.0 * pi;
        double magic = Math.Sin(radLat);
        magic = 1 - ee * magic * magic;
        double sqrtMagic = Math.Sqrt(magic);
        dLat = (dLat * 180.0) / ((a * (1 - ee)) / (magic * sqrtMagic) * pi);
        dLon = (dLon * 180.0) / (a / sqrtMagic * Math.Cos(radLat) * pi);
        double mgLat = lat + dLat;
        double mgLon = lon + dLon;
        return new Vector2D(mgLat, mgLon);
    }

    /// <summary>
    /// 火星坐标系 (GCJ-02) to GPS84
    /// </summary>
    /// <param name="lat"></param>
    /// <param name="lon"></param>
    /// <returns></returns>
    public Vector2D gcj02_To_Gps84(double lat, double lon)
    {
        Vector2D gps = gps84_To_Gcj02(lat, lon);
        double lontitude = lon * 2 - gps.y;
        double latitude = lat * 2 - gps.x;
        return new Vector2D(latitude, lontitude);
    }


    /// <summary>
    /// 将 GCJ-02 坐标转换成 BD-09 坐标
    /// </summary>
    /// <param name="lat"></param>
    /// <param name="lon"></param>
    /// <returns></returns>
    public Vector2D gcj02_To_Bd09(double lat, double lon)
    {
        double x = lon, y = lat;
        double z = Math.Sqrt(x * x + y * y) + 0.00002 * Math.Sin(y * x_pi);
        double theta = Math.Atan2(y, x) + 0.000003 * Math.Cos(x * x_pi);
        double tempLon = z * Math.Cos(theta) + 0.0065;
        double tempLat = z * Math.Sin(theta) + 0.006;
        return new Vector2D(tempLat, tempLon);
    }

    /// <summary>
    /// 百度BD-09转火星坐标系GCJ-02
    /// </summary>
    /// <param name="lat"></param>
    /// <param name="lon"></param>
    /// <returns></returns>
    public Vector2D bd09_To_Gcj02(double lat, double lon)
    {
        double x = lon - 0.0065, y = lat - 0.006;
        double z = Math.Sqrt(x * x + y * y) - 0.00002 * Math.Sin(y * x_pi);
        double theta = Math.Atan2(y, x) - 0.000003 * Math.Cos(x * x_pi);
        double tempLon = z * Math.Cos(theta);
        double tempLat = z * Math.Sin(theta);
        double[] gps = { tempLat, tempLon };
        return new Vector2D(tempLat, tempLon);
    }

    /// <summary>
    ///  将GPS84转为百度BD-09
    /// </summary>
    /// <param name="lat"></param>
    /// <param name="lon"></param>
    /// <returns></returns>
    public Vector2D gps84_To_bd09(double lat, double lon)
    {
        Vector2D gcj02 = gps84_To_Gcj02(lat, lon);
        Vector2D bd09 = gcj02_To_Bd09(gcj02.x, gcj02.y);
        return bd09;
    }

    /// <summary>
    /// 百度BD09转GPS84
    /// </summary>
    /// <param name="lat"></param>
    /// <param name="lon"></param>
    /// <returns></returns>
    public Vector2D bd09_To_gps84(double lat, double lon)
    {
        Vector2D gcj02 = bd09_To_Gcj02(lat, lon);
        Vector2D gps84 = gcj02_To_Gps84(gcj02.x, gcj02.y);
        gps84.x = Math.Round(gps84.x, 6);
        gps84.y = Math.Round(gps84.y, 6);
        return gps84;
    }

    #endregion
}
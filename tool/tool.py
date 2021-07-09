# latitude:纬度 longitude:经度 altitude:海拔
import math
import numpy


###############################################################################################
# 工具类，包含：
#        将经纬度转换成大地坐标系XYZ
#        将大地坐标系转换成经纬度LLA
###############################################################################################
def LLA_to_XYZ(latitude, longitude, altitude):
    # 经纬度的余弦值
    cosLat = math.cos(latitude * math.pi / 180)
    sinLat = math.sin(latitude * math.pi / 180)
    cosLon = math.cos(longitude * math.pi / 180)
    sinLon = math.sin(longitude * math.pi / 180)

    # WGS84坐标系的参数
    rad = 6378137.0  # 地球赤道平均半径（椭球长半轴：a）
    f = 1.0 / 298.257224  # WGS84椭球扁率 :f = (a-b)/a
    C = 1.0 / math.sqrt(cosLat * cosLat + (1 - f) * (1 - f) * sinLat * sinLat)
    S = (1 - f) * (1 - f) * C
    h = altitude

    # 计算XYZ坐标
    X = (rad * C + h) * cosLat * cosLon
    Y = (rad * C + h) * cosLat * sinLon
    Z = (rad * S + h) * sinLat

    return numpy.array([X, Y, Z])


def XYZ_to_LLA(X, Y, Z):
    # WGS84坐标系的参数
    a = 6378137.0  # 椭球长半轴
    b = 6356752.314245  # 椭球短半轴
    ea = numpy.sqrt((a ** 2 - b ** 2) / a ** 2)
    eb = numpy.sqrt((a ** 2 - b ** 2) / b ** 2)
    p = numpy.sqrt(X ** 2 + Y ** 2)
    theta = numpy.arctan2(Z * a, p * b)
    # 计算经纬度及海拔
    longitude = numpy.arctan2(Y, X)
    latitude = numpy.arctan2(Z + eb ** 2 * b * numpy.sin(theta) ** 3, p - ea ** 2 * a * numpy.cos(theta) ** 3)
    N = a / numpy.sqrt(1 - ea ** 2 * numpy.sin(latitude) ** 2)
    altitude = p / numpy.cos(latitude) - N
    res = list(numpy.array([numpy.degrees(latitude), numpy.degrees(longitude), altitude]))
    return res

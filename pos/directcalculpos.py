import numpy as np

#########################################################################################
# 通过方程转化直接求出最后坐标XYZ和已知参数的关系，直接将参数带入求解
#########################################################################################
inf = [[-368461.739, 26534822.568, -517664.322, 21966984.2427, -0.000104647296],
       [10002180.758, 12040222.131, 21796269.831, 23447022.1136, -0.000308443058],
       [-7036480.928, 22592611.906, 11809485.040, 20154521.4618, -0.000038172460],
       [8330122.410, 23062955.196, 10138101.718, 22129309.3677, -0.0002393560]]

inf1 = [[-368461.739, 26534822.568, -517664.322],
        [10002180.758, 12040222.131, 21796269.831],
        [-7036480.928, 22592611.906, 11809485.040],
        [8330122.41, 23062955.196, 10138101.718]]

p = [21966984.2427, 23447022.1136, 20154521.4618, 22129309.3677]

t = [-0.000104647296, -0.000308443058, -0.000038172460, -0.000239356]
K = []
r = []
C = []
x = []
y = []
z = []
pos = []
pos1 = [-2279828.8292, 5004706.5483, 3219777.4684]


def calKi():
    i = 0
    j = len(inf1)
    while i < j:
        r = round((inf1[i][0] ** 2), 5) + round((inf1[i][1] ** 2), 5) + round((inf1[i][2] ** 2), 5)
        # print((inf1[i][0] ** 2))
        # print((inf1[i][1] ** 2))
        # print((inf1[i][2] ** 2))
        K.append(r)
        i = i + 1
    return K


def calr():
    global r
    i = 0
    j = len(p)
    res = 0
    while i < j:
        res = p[i] + (3 * 10 ** 8) * t[i]
        r.append(res)
        i = i + 1
    return r


def calC():
    global C
    i = 1
    j = len(p)
    while i < j:
        re = r[i] ** 2 - r[0] ** 2
        res = K[i] - K[0] - re
        C.append(res)
        i = i + 1
    return C


def calxyz():
    global x
    global y
    global z
    i = 1
    j = len(p)
    while i < j:
        x.append(round(inf1[i][0] - inf1[0][0], 3))
        y.append(round(inf1[i][1] - inf1[0][1], 3))
        z.append(round(inf1[i][2] - inf1[0][2], 3))
        i = i + 1


def XYZ_to_LLA(X, Y, Z):
    # WGS84坐标系的参数
    a = 6378137.0  # 椭球长半轴
    b = 6356752.314245  # 椭球短半轴
    ea = np.sqrt((a ** 2 - b ** 2) / a ** 2)
    eb = np.sqrt((a ** 2 - b ** 2) / b ** 2)
    p = np.sqrt(X ** 2 + Y ** 2)
    theta = np.arctan2(Z * a, p * b)
    # 计算经纬度及海拔
    longitude = np.arctan2(Y, X)
    latitude = np.arctan2(Z + eb ** 2 * b * np.sin(theta) ** 3, p - ea ** 2 * a * np.cos(theta) ** 3)
    N = a / np.sqrt(1 - ea ** 2 * np.sin(latitude) ** 2)
    altitude = p / np.cos(latitude) - N
    res = np.array([np.degrees(latitude), np.degrees(longitude), altitude])
    return list(res)


def calpos():
    global pos
    dd = x[0] * y[1] * z[2] + y[0] * z[1] * x[2] + x[1] * y[2] * z[0] - z[0] * y[1] * x[2] - y[0] * x[1] * z[2] - z[1] * \
         y[2] * x[0]
    D = 2 * dd
    print("D:", D)
    D1 = C[0] * (y[1] * z[2] - y[2] * z[1]) + C[1] * (y[2] * z[0] - y[0] * z[2]) + C[2] * (y[0] * z[1] - y[1] * z[0])
    D2 = C[0] * (x[2] * z[1] - x[1] * z[2]) + C[1] * (x[0] * z[2] - x[2] * z[0]) + C[2] * (x[1] * z[0] - x[0] * z[1])
    D3 = C[0] * (x[1] * y[2] - x[2] * y[1]) + C[1] * (x[2] * y[0] - x[0] * y[2]) + C[2] * (x[0] * y[1] - x[1] * y[0])
    print("D1,D2,D3:", D1, D2, D3)
    X1 = D1 / D
    Y1 = D2 / D
    Z1 = D3 / D
    pos.append(X1)
    pos.append(Y1)
    pos.append(Z1)
    detpos = np.array(pos) - np.array(pos1)
    print("目标位置的XYZ：", pos1)
    print("目标位置的经纬度：", XYZ_to_LLA(pos1[0], pos1[1], pos1[2]))
    print("直接计算的XYZ坐标：", pos)
    print("直接计算的经纬度：", XYZ_to_LLA(pos[0], pos[1], pos[2]))
    print("计算位置和目标位置的距离差：", detpos)


calKi()
print("K:", K)
calr()
print("r:", r)
calC()
print("C:", C)
calxyz()
print("x:", x)
print("y:", y)
print("z:", z)
calpos()

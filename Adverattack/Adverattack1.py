import math
import torch
import numpy as np
import csv
from scipy.optimize import fsolve

f = open('../dataset/mydataset.csv', 'r')
reader = csv.reader(f)  # 创建一个与该文件相关的阅读器
result = list(reader)
inf = []
inf1 = []
p = []
t = []
attpos = [-2280702.2803424736, 5004308.566575137, 3219777.468027828]
# 初始估计接收机位置
pos0 = [0, 0, 0]
x1 = []
y1 = []
z1 = []
Ts = 0


# 从文件中获取数据
def getdata(i):
    global inf
    k = i
    j = i + 4
    while k < j:
        a = [float(result[k][2]), float(result[k][3]), float(result[k][4]), float(result[k][5]), float(result[k][6])]
        inf.append(a)
        k += 1
    return inf


# 将获取的数据存入到相应的数组中
def processdata():
    global inf1
    global p
    global t
    for i in range(4):
        r = []
        r.append(inf[i][0])
        r.append(inf[i][1])
        r.append(inf[i][2])
        inf1.append(r)
        p.append(inf[i][3])
        t.append(inf[i][4])


# 获取两点之间的距离
def get_distance1(pos1, pos2):
    v1 = np.array([pos1[0], pos1[1], pos1[2]])
    v2 = np.array([pos2[0], pos2[1], pos2[2]])
    distance = np.linalg.norm(v1 - v2)
    return distance


# 获取修复钟差的伪距
def get_distance2(pos1, pos2, clk):
    v1 = np.array([pos1[0], pos1[1], pos1[2]])
    v2 = np.array([pos2[0], pos2[1], pos2[2]])
    distance = np.linalg.norm(v1 - v2)
    res = distance - (3 * 10 ** 8) * clk
    return res


# 获取估计伪距值
def getp():
    res = []
    i = 0
    while i < len(inf1):
        res.append(get_distance2(pos0, inf1[i], t[i]))
        i = i + 1
    return res


# 获取伪距差值
def getdetp(p1):
    res = []
    k = len(p)
    i = 0
    while i < k:
        res.append(p[i] - p1[i])
        i = i + 1
    return res


# 获得每一个r
def getdis():
    res = []
    k = len(p)
    i = 0
    while i < k:
        dis = get_distance1(inf[i], pos0)
        res.append(dis)
        i = i + 1
    return res


# 获取观测矩阵
def getmatH(info, pos, r):
    res = []
    k = len(info)
    i = 0
    while i < k:
        l = []
        l.append((pos[0] - info[i][0]) / r[i])
        l.append((pos[1] - info[i][1]) / r[i])
        l.append((pos[2] - info[i][2]) / r[i])
        i = i + 1
        res.append(l)
    return res


def calTs(p1, H):
    global pos0
    detp1 = getdetp(p1)
    Ho = H
    # print(list(Ho))
    H1 = np.array(Ho)
    # print("观测矩阵:\n", H1)
    H2 = np.transpose(H1)
    H3 = np.dot(H2, H1)
    # print(H3)
    H4 = np.linalg.pinv(H3)
    H5 = np.dot(H1, H4)
    H6 = np.dot(H5, H2)
    # print(H6)
    In = np.eye(len(p))
    S = In - H6
    R = np.dot(S, detp1)
    Rt = np.transpose(R)
    Ts = np.dot(Rt, R)
    # print("Ts:", Ts)
    return Ts


# 计算接收机位置
def calresult():
    global pos0
    global Ts
    for j in range(10):
        # print(pos0)
        p1 = getp()
        # print("估计位置到各卫星的伪距值为：\n", p1)
        detp = getdetp(p1)
        # print("估计伪距和实际伪距的差：\n", detp)
        r = getdis()
        # print("获取的r为：\n", r)
        H = getmatH(inf1, pos0, r)
        # print("获得的观测矩阵H为:\n", H)
        Ts = calTs(p1, H)
        H1 = np.array(H)
        # print("观测矩阵:\n", H1)
        H2 = np.transpose(H1)
        # print(H2)
        H3 = np.dot(H2, H1)
        # print(H3)
        H4 = np.linalg.pinv(H3)
        # print(H4)
        H5 = np.dot(H4, H2)
        detx = np.dot(H5, detp)
        pos0 = pos0 + detx


# 将XYZ坐标转换成经纬度
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
    return np.array([np.degrees(latitude), np.degrees(longitude), altitude])


# 获取梯度矩阵
def getgrad1():
    global p
    last = []
    for j in range(4):
        # print(j)
        for i in range(10):
            p[j] = p[j] + i / 10
            last = pos0
            calresult()
            # print("pos0", pos0)
            # print("last", last)
        if pos0[0] < last[0]:
            x1.append(1)
        else:
            x1.append(-1)
        if pos0[1] < last[1]:
            y1.append(1)
        else:
            y1.append(-1)
        if pos0[2] < last[2]:
            z1.append(1)
        else:
            z1.append(-1)
    # print(x1)
    # print(y1)
    # print(z1)


# 获取梯度符号
def getgrad(x, y, z):
    n1 = x - attpos[0]
    n2 = y - attpos[1]
    n3 = z - attpos[2]
    if n1 > 0:
        gx = np.array(x1)
    else:
        gx = -np.array(x1)
    if n2 > 0:
        gy = np.array(y1)
    else:
        gy = -np.array(y1)
    if n3 > 0:
        gz = np.array(z1)
    else:
        gz = -np.array(z1)
    g = np.array(gx) + np.array(gy) + np.array(gz)
    # t = torch.tensor(g)
    # sign = t.sign()
    return g


# 对抗攻击方法
def adattack():
    global inf1
    global p
    global t
    global Ts
    dis = 1000000
    while dis > 5:
        calresult()
        r = pos0
        print("攻击位置：", attpos)
        print("解算位置：", r)
        print("经纬度：", XYZ_to_LLA(r[0], r[1], r[2]))
        dis = get_distance1(r, attpos)
        print("距离差：", dis)
        print("TS:", Ts)
        sign = getgrad(r[0], r[1], r[2])
        # print(sign)
        i = 1
        while i < 4:
            p[i] = p[i] + 1 * sign[i]
            i = i + 1


getdata(21)
processdata()
getgrad1()
adattack()

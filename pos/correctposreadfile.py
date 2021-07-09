import csv
import tool.tool as tool
import numpy as np

######################################################################################################
# 测试版伪距迭代定位，数据从文件中读取，使用的是九峰站的数据
######################################################################################################
f = open('../dataset/mydataset.csv', 'r')
reader = csv.reader(f)  # 创建一个与该文件相关的阅读器
result = list(reader)
inf = []
inf1 = []
p = []
t = []
# 初始估计接收机位置
pos0 = [0, 0, 0]
pos1 = [-2279829.1069, 5004709.2387, 3219779.0559]


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
def get_distance1(pos3, pos2):
    v1 = np.array([pos3[0], pos3[1], pos3[2]])
    v2 = np.array([pos2[0], pos2[1], pos2[2]])
    distance = np.linalg.norm(v1 - v2)
    return distance


# 获取修复钟差的伪距
def get_distance2(pos3, pos2, clk):
    v1 = np.array([pos3[0], pos3[1], pos3[2]])
    v2 = np.array([pos2[0], pos2[1], pos2[2]])
    distance = np.linalg.norm(v1 - v2)
    res = distance - (3 * 10 ** 8) * clk
    return res


# 获取估计位置到卫星的估计伪距值
def getp():
    res = []
    i = 0
    while i < len(inf1):
        res.append(get_distance2(pos0, inf1[i], t[i]))
        i = i + 1
    return res


# 获取估计伪距与观测伪距之间的差值
def getdetp(p1):
    res = []
    k = len(p)
    i = 0
    while i < k:
        res.append(p[i] - p1[i])
        i = i + 1
    return res


# 获得接收机估计位置到每一颗卫星的几何距离
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
    c = 3 * 10 ** 8
    while i < k:
        l = []
        l.append((pos[0] - info[i][0]) / r[i])
        l.append((pos[1] - info[i][1]) / r[i])
        l.append((pos[2] - info[i][2]) / r[i])
        i = i + 1
        res.append(l)
    return res


# 计算接收机的估计位置
def calresult():
    global pos0
    for j in range(10):
        p1 = getp()
        # print("估计位置到各卫星的伪距值为：\n", p1)
        detp = getdetp(p1)
        # print("估计伪距和实际伪距的差：\n", detp)
        r = getdis()
        # print("获取的r为：\n", r)
        H = getmatH(inf1, pos0, r)
        # print("获得的观测矩阵H为:\n", H)
        H1 = np.array(H)
        # print(H1)
        H2 = np.transpose(H1)
        # print(H2)
        H3 = np.dot(H2, H1)
        # print(H3)
        H4 = np.linalg.pinv(H3)
        # print(H4)
        H5 = np.dot(H4, H2)
        # print(H5)
        det = np.dot(H5, detp)
        # print(det)
        pos0 = pos0 + det
        print("接收机估计位置为：；", pos0)
    # print(pos1 - pos0)


getdata(17)
processdata()
calresult()
print("接收机经纬度位置", tool.XYZ_to_LLA(pos0[0], pos0[1], pos0[2]))

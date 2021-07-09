import numpy as np
import tool.tool as tool

######################################################################################################
# 测试版RAIM,使用的是BDS拉萨站的数据
######################################################################################################
inf = [[-14760528.492, 39470117.279, -256073.994, 36936614.601, 0.000566140763],
       [-10925771.656, 23917731.411, 33139431.442, 36762486.288, 0.000135112551],
       [-4940290.027, 41376945.710, 5236336.815, 35957376.805, 0.000852689119],
       [6294727.557, 31138188.359, 28056129.289, 36542318.068, -0.000857045235]]

inf1 = [[-14760528.492, 39470117.279, -256073.994],
        [-10925771.656, 23917731.411, 33139431.442],
        [-4940290.027, 41376945.710, 5236336.815],
        [6294727.557, 31138188.359, 28056129.289]]

p = [36936614.601, 36762486.288, 35957376.805, 36542318.068]

t = [0.000566140763, 0.000135112551, 0.000852689119, -0.000857045235]
# 初始估计接收机位置
pos0 = [0, 0, 0]

pos1 = [-106943.5, 5549296.14, 3139212.6]
x1 = []
y1 = []
z1 = []


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


# 获取观测矩阵
def getmatH2(info, pos, r):
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


def calTs(p1, H):
    global pos0
    detp1 = getdetp(p1)
    # print(list(Ho))
    H1 = np.array(H)
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
    print("Ts:", Ts)


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
        calTs(p1, H)
        # print(H5)
        det = np.dot(H5, detp)
        # print(det)
        pos0 = pos0 + det
        print("接收机估计位置为：；", pos0)
    print(pos1 - pos0)


# 获取攻击使用的梯度矩阵
def getgrad():
    global p
    last = []
    for j in range(4):
        print(j)
        for i in range(10):
            p[j] = p[j] + i / 10
            last = pos0
            calresult()
            print("pos0", pos0)
            print("last", last)
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
    print(x1)
    print(y1)
    print(z1)


calresult()
print(tool.XYZ_to_LLA(pos0[0], pos0[1], pos0[2]))
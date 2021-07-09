import tool.tool as tool

pos = [30.51556540968659, 114.50101889963045, 71.32352626137435]
print(list(tool.LLA_to_XYZ(pos[0], pos[1], pos[2])))

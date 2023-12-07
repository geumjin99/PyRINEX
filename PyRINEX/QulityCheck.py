import json
from PyRINEX.reader import *
from PyRINEX.DataManagement import *
import numpy as np
import math
from numpy import *
import csv
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from dateutil import parser
import time
from datetime import datetime
import multiprocessing
import threading
import concurrent.futures
import seaborn as sns
def xyz2BLH(x, y, z, rad=True):#output is lon lat. This function is from https://zhuanlan.zhihu.com/p/358210424
    a = 6378137.0000
    b = 6356752.3141
    e2 = 1 - (b / a)**2
    p = np.sqrt(x**2+y**2)
    L = np.arctan2(y, x)
    def rad2degree(rad):
        return rad * 180 / np.pi
    def ite(z,p):
        def cal_N(B): return a/np.sqrt(1-e2*np.sin(B)**2)
        def cal_H(N, B): return p/np.cos(B)-N
        def cal_B(N, H): return np.arctan(z/((1 - e2*N/(N+H))*p))
        B = cal_B(1, 0)
        N = cal_N(B)
        H0, H = 1e9, cal_H(N, B)
        while np.abs(H - H0) > 0.1:
            B = cal_B(N, H)
            N = cal_N(B)
            H0, H = H, cal_H(N, B)
        return H, B
    H, B = np.vectorize(ite)(z,p)
    if rad:
        return [L, B * 1, H * 1]
    else:
        return [rad2degree(L), rad2degree(B), H * 1]
def utf2jd(utf): #This function is from https://zhuanlan.zhihu.com/p/435129629
    year = utf[0]
    if len(year) == 2:
        year = "20" + year
    year = int(year)
    month = utf[1]
    day = int(utf[2]) + (int(utf[3]) * 3600 + int(utf[4]) * 60 + float(utf[5])) / 86400
    y = int(year) + 4800
    m = int(month)
    if m <= 2:
        m = m + 12
        y = y - 1
    e = math.floor(30.6 * (m+1))
    a = math.floor(y/100)
    if (year < 1582) or (year == 1582 and month < 10) or (year == 1582 and month == 10 and day < 15):
        b = -38
    else:
        b = math.floor((a / 4) - a)
    c = math.floor(365.25 * y)
    jd = b + c + e + day - 32167.5
    mjd = jd - 2400000.5
    return mjd
def utf2gps(utf):
    mjd = utf2jd(utf)
    e = mjd - 44244
    week = math.floor(e/7)
    e = e - week*7
    return [week, round(e*86400)]
def t_kCal(gps_seconds, toe):
    t_k = gps_seconds - toe
    if t_k > 302400:
        t_k -= 604800
    if t_k < -302400:
        t_k += 604800
    if -7200 <= t_k and t_k <= 7200:
        t_k = t_k
    else:
        t_k = "not valid"
    return t_k
def SatelliteSignalPlot(opath):
    odata = json.loads(observations(opath))
    obss = odata
    id = []
    TIME = []
    for k,v in obss.items():
        for ks,vs in obss[k].items():
            if type(obss[k][ks]) == dict:
                id.append(ks)
                TIME.append(k.split()[2] + " " + ":".join(k.split()[3:6]))
    ax = plt.subplot()
    hours = mdates.DateFormatter('%d %H:%M')
    ax.xaxis.set_major_formatter(hours)
    ax.xaxis.set_major_formatter(hours)
    ax.set_title(os.path.basename(opath)[:-4] + " Signal Plot")
    time = [parser.parse(x) for x in TIME]
    plt.scatter(time, id, s=2, marker=".", alpha=0.75, c="g")
    plt.grid(linestyle = '--', linewidth=0.5)
    plt.ylabel("SATELLITE NO")
    plt.gcf().set_size_inches(15, 7)
    plt.yticks(fontsize=5)
    filename = opath[:-4] + "SignalPlot.jpg"
    plt.savefig(filename, dpi=360)
    plt.cla()
    plt.clf()
def time_serises(odata):
    obss = odata
    list = []
    for k,v in obss.items():
        list.append(k)
    return list
def plot(filename, title, gpsserises, epochs, dataset, type, column, limit=100, y_label=None):
    TIME = [epoch.split()[2] + " " + ":".join(epoch.split()[3:6]) for epoch in epochs]
    hours = mdates.DateFormatter('%d %H:%M')
    time = [parser.parse(x) for x in TIME]

    ax = plt.subplot()
    hours = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(hours)
    ax.set_title(title)

    element = {}
    gps_prn = gpsserises
    number_of_plots = len(gps_prn)

    # 使用seaborn中的颜色调色板
    colors = sns.color_palette("viridis", number_of_plots)
    ax.set_prop_cycle('color', colors)

    for n in range(len(gps_prn)):
        element[gps_prn[n]] = []

    for n in range(len(dataset)):
        for m in range(len(gps_prn)):
            if dataset[n][m][type] != 0 and abs(dataset[n][m][type]) < limit:
                element[gps_prn[m]].append([time[n], dataset[n][m][type]])

    for n in range(len(gps_prn)):
        data = element[gps_prn[n]]
        x = [entry[0] for entry in data]
        y = [entry[1] for entry in data]
        ax.scatter(x, y, s=1, label=gps_prn[n])

    plt.legend(bbox_to_anchor=(1.05, 0), loc=3, borderaxespad=0, ncol=column)

    if y_label:
        ax.set_ylabel(y_label)

    plt.savefig(filename, dpi=360, bbox_inches='tight')
    plt.cla()
def gps_pos(opath):
    odata = json.loads(observations(opath))
    if opath[-1] == "o":
        npath = opath[:-1] + "n"
    else:
        npath = opath[:-1] + "N"
    ndata = json.loads(navigations(npath))
    pos_list = []
    a = 6378137
    f = 1/298.257223563
    Omega_e_Dot = 7.2921151467e-5
    mu = 3.986005e14
    c = 2.99792458e8
    epoch_index = odata
    postions = []
    for k,v in epoch_index.items():
        postions_per_epoch = {}
        postions_per_epoch["epoch"] = k
        gps = []
        sondic = epoch_index[k]
        for ks, vs in sondic.items():
            utf = k.split()
            gps_seconds = utf2gps(utf)[1]
            if "G" in ks:
                gps.append(ks)

        for x in range(len(gps)):
            index = []
            t_kResult = []
            for m in range(len(ndata[gps[x]])):
                toe = float(ndata[gps[x]][m]["Toe Time of Ephemeris"].replace("D", "E"))
                try:
                    if t_kCal(gps_seconds, toe) != "not valid":
                        index.append(m)
                        t_kResult.append(abs(t_kCal(gps_seconds, toe)))
                except ValueError:
                    index.append(m)
                    t_kResult.append(abs(t_kCal(gps_seconds, toe)))
            dic = ndata[gps[x]][index[t_kResult.index(min(t_kResult))]]
            toe = float(dic["Toe Time of Ephemeris"].replace("D", "E"))
            A_sqrt = float(dic["sqrt(A)"].replace("D", "E"))
            e = float(dic["e Eccentricity"].replace("D", "E"))
            i_0 = float(dic["i0"].replace("D", "E"))
            Omega_0 = float(dic["OMEGA0"].replace("D", "E"))
            omega = float(dic["omega"].replace("D", "E"))
            M_0 = float(dic["M0"].replace("D", "E"))
            Delta_n = float(dic["Delta n"].replace("D", "E"))
            i_Dot = float(dic["IDOT"].replace("D", "E"))
            Omega_Dot = float(dic["OMEGA DOT"].replace("D", "E"))
            Cuc = float(dic["Cuc"].replace("D", "E"))
            Cus = float(dic["Cus"].replace("D", "E"))
            Crc = float(dic["Crc"].replace("D", "E"))
            Crs = float(dic["Crs"].replace("D", "E"))
            Cic = float(dic["Cic"].replace("D", "E"))
            Cis = float(dic["CIS"].replace("D", "E"))
            A = A_sqrt ** 2
            t_k = t_kCal(gps_seconds, toe)
            n_0 = math.sqrt(mu/A**3)
            n = n_0 + Delta_n
            M_k = M_0 + n*t_k
            if M_k < 0:
                M_k += 2 * math.pi
            if M_k > 2 * math.pi:
                M_k -= 2 * math.pi
            E_old = M_k
            E_new = M_k + e * math.sin(E_old)
            i = 1
            while abs(E_new - E_old) > 1e-8:
                E_old = E_new
                E_new = M_k + e * math.sin(E_old)
                i += 1
                if (i > 10):
                    break
            E_k = E_new
            cosNu_k = (math.cos(E_k) - e) / (1 - e * math.cos(E_k))
            sinNu_k = (math.sqrt(1 - e ** 2) * math.sin(E_k)) / (1 - e * math.cos(E_k))
            if cosNu_k == 0:
                if sinNu_k > 0:
                    Nu_k = math.pi / 2
                else:
                    Nu_k = -math.pi / 2
            else:
                Nu_k = math.atan(sinNu_k / cosNu_k)
            if cosNu_k < 0:
                if sinNu_k >= 0:
                    Nu_k += math.pi
                else:
                    Nu_k -= math.pi
            Phi_k = Nu_k + omega
            delta_u_k = Cus * math.sin(2 * Phi_k) + Cuc * math.cos(2 * Phi_k)
            delta_r_k = Crs * math.sin(2 * Phi_k) + Crc * math.cos(2 * Phi_k)
            delta_i_k = Cis * math.sin(2 * Phi_k) + Cic * math.cos(2 * Phi_k)
            u_k = Phi_k + delta_u_k
            r_k = A * (1 - e * math.cos(E_k)) + delta_r_k
            i_k = i_0 + i_Dot * t_k + delta_i_k
            x_p_k = r_k * math.cos(u_k)
            y_p_k = r_k * math.sin(u_k)
            Omega_k = Omega_0 + (Omega_Dot - Omega_e_Dot) * t_k - Omega_e_Dot * toe
            x_k = x_p_k * math.cos(Omega_k) - y_p_k * math.cos(i_k) * math.sin(Omega_k)
            y_k = x_p_k * math.sin(Omega_k) + y_p_k * math.cos(i_k) * math.cos(Omega_k)
            z_k = y_p_k * math.sin(i_k)
            postions_per_epoch[gps[x]] = [x_k, y_k, z_k]
        postions.append(postions_per_epoch)
    return postions
def azi_ele(opath):
    header = json.loads(oheader(opath))
    odata = json.loads(observations(opath))
    if opath[-1] == "o":
        npath = opath[:-1] + "n"
    else:
        npath = opath[:-1] + "N"
    ndata = json.loads(navigations(npath))
    prns = header["PRNS"]
    epoch_index = odata
    xyz = header["APPROX POSITION XYZ"]
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    p = math.sqrt(x ** 2 + y ** 2)
    R = math.sqrt(x ** 2 + y ** 2 + z ** 2)
    e = [[-y / p, x / p, 0],
        [-x * z / (p * R), -y * z / (p * R), p / R],
        [x / R, y / R, z / R]]
    postion = gps_pos(opath)
    azi_ele_list = []
    for n in range(len(postion)):
        dic = {}
        dic["epoch"] = postion[n]["epoch"]
        gps = []
        for k, v in postion[n].items():
            if "G" in k:
               gps.append(k)
        for j in range(len(gps)):
            pos = postion[n][gps[j]]
            d = [0, 0, 0]
            for m in range(3):
                d[m] = 0
                for o in range(3):
                    d[m] = d[m] + (pos[o] - xyz[o])*(e[m][o])
            s = d[2]/(math.sqrt(d[0]**2 + d[1]**2 + d[2]**2))
            if s == 1:
                ele = 0.5 * math.pi
            else:
                ele = math.atan(s / math.sqrt(1 - s ** 2))
            if d[1] == 0 and d[0] > 0:
                azi = 0.5 * math.pi
            elif d[1] == 0 and d[0] < 0:
                azi = 1.5 * math.pi
            else:
                azi = math.atan(d[0] / d[1])
                if d[1] < 0:
                    azi = azi + math.pi
                elif d[1] > 0 and d[0] < 0:
                    azi = azi + 2 * math.pi
            dic[gps[j]] = [(azi), (ele)]
        azi_ele_list.append(dic)
#    ax = plt.subplot(111, polar=True)
#    ax.set_theta_direction(-1)
#    ax.set_theta_zero_location('N')
#    ax.set_rlim(90, 0)
#    ax.set_title(os.path.basename(opath)[:-4]+"SkyPlot")
#    element = {}
#    gpss = []
#    for n in range(len(prns)):
#        if "G" in prns[n]:
#            gpss.append(prns[n])
#    colormap = plt.cm.nipy_spectral
#    colors = [colormap(i) for i in np.linspace(0, 1, len(gpss))]
#    ax.set_prop_cycle('color', colors)
#    for n in range(len(gpss)):
#        element[gpss[n]] = []
#    for n in range(len(azi_ele_list)):
#        for k, v in azi_ele_list[n].items():
#            for keys, values in element.items():
#                if k == keys:
#                    element[keys].append(azi_ele_list[n][k])
#    for n in range(len(gpss)):
#        data = element[gpss[n]]
#        rho = []
#        theta = []
#        for m in range(len(data)):
#            rho.append(data[m][0])
#            theta.append(math.degrees(data[m][1]))
#        ax.scatter(rho, theta, s=10, label=gpss[n])
##
#    plt.legend(bbox_to_anchor=(0, 1.05), borderaxespad=3)
#    plt.savefig(opath[:-4]+"skyplot.jpg", dpi=360, bbox_inches='tight')
#    plt.cla()
#    plt.clf()
#    fieldnames = gpss
#    fieldnames.insert(0, "epoch")
#    row = []
#    for n in range(len(azi_ele_list)):
#        element = {}
#        element["epoch"] = azi_ele_list[n]["epoch"]
#        for k, v in azi_ele_list[n].items():
#            element[k] = azi_ele_list[n][k]
#        row.append(element)
    ax = plt.subplot(111, polar=True)
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N')
    ax.set_rlim(90, 0)
    ax.set_title(os.path.basename(opath)[:-4] + "SkyPlot")

    element = {}
    gpss = [prn for prn in prns if "G" in prn]

    # 使用 'viridis' 调色板，避免红色和绿色
    colormap = plt.cm.viridis
    colors = [colormap(i) for i in np.linspace(0, 1, len(gpss))]
    ax.set_prop_cycle('color', colors)

    for n in range(len(gpss)):
        element[gpss[n]] = []

    for n in range(len(azi_ele_list)):
        for k, v in azi_ele_list[n].items():
            for keys, values in element.items():
                if k == keys:
                    element[keys].append(azi_ele_list[n][k])

    for n in range(len(gpss)):
        data = element[gpss[n]]
        rho = []
        theta = []
        for m in range(len(data)):
            rho.append(data[m][0])
            theta.append(math.degrees(data[m][1]))
        ax.scatter(rho, theta, s=10, label=gpss[n])

    plt.legend(bbox_to_anchor=(0, 1.05), borderaxespad=3)
    plt.savefig(opath[:-4] + "skyplot.png", dpi=360, bbox_inches='tight')
    plt.cla()
    plt.clf()

    # 将数据保存到CSV文件
    fieldnames = gpss
    fieldnames.insert(0, "epoch")
    row = []
    for n in range(len(azi_ele_list)):
        element = {"epoch": azi_ele_list[n]["epoch"]}
        for k, v in azi_ele_list[n].items():
            element[k] = azi_ele_list[n][k]
        row.append(element)
    with open(opath[:-4] +"aziele.csv", "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames = fieldnames)
        writer.writeheader()
        writer.writerows(row)
    return azi_ele_list
def LandC(opath):
    odata = json.loads(observations(opath))
    prns = json.loads(oheader(opath))["PRNS"]
    amount = len(prns)
    time_index = time_serises(odata)
    epoch_index = odata
    size = len(epoch_index)
    lc = np.zeros(shape=(size, amount, 4), dtype=float)
    index = -1
    for k, v in epoch_index.items():
        sondic = epoch_index[k]
        index +=1
        for ks, vs in sondic.items():
            if "G" in ks:
                sondic1 = sondic[ks]
                obs = []
                for kss, vss in sondic1.items():
                    if "L1" in kss and sondic1[kss] != "         0.000  ":
                        obs.append(float(sondic1[kss][0:14]))
                if len(obs) != 1:
                    if len(obs) == 0:
                        obs.append(0)
                    if len(obs) > 1:
                        while len(obs) > 1:
                            obs.pop()
                for kss, vss in sondic1.items():
                    if "L2" in kss and sondic1[kss] != "         0.000  ":
                        obs.append(float(sondic1[kss][0:14]))
                if len(obs) != 2:
                    if len(obs) == 1:
                        obs.append(0)
                    if len(obs) > 2:
                        while len(obs) > 2:
                            obs.pop()
                for kss, vss in sondic1.items():
                    if "C1" in kss and sondic1[kss] != "         0.000  ":
                        obs.append(float(sondic1[kss][0:14]))
                    if "P1" in kss and sondic1[kss] != "         0.000  ":
                        obs.append(float(sondic1[kss][0:14]))
                if len(obs) != 3:
                    if len(obs) == 2:
                        obs.append(0)
                    if len(obs) > 3:
                        while len(obs) > 3:
                            obs.pop()
                for kss, vss in sondic1.items():
                    if "C2" in kss and sondic1[kss] != "         0.000  ":
                        obs.append(float(sondic1[kss][0:14]))
                    if "P2" in kss and sondic1[kss] != "         0.000  ":
                        obs.append(float(sondic1[kss][0:14]))
                if len(obs) != 4:
                    if len(obs) == 3:
                        obs.append(0)
                    if len(obs) > 4:
                        while len(obs) > 4:
                            obs.pop()
                lc[time_index.index(k)][prns.index(ks)][0] = obs[0]
                lc[time_index.index(k)][prns.index(ks)][1] = obs[1]
                lc[time_index.index(k)][prns.index(ks)][2] = obs[2]
                lc[time_index.index(k)][prns.index(ks)][3] = obs[3]
            if "R" in ks:
                sondic1 = sondic[ks]
                obs = []
                for kss, vss in sondic1.items():
                    if "L1" in kss and sondic1[kss] != "         0.000  ":
                        obs.append(float(sondic1[kss][0:14]))
                if len(obs) != 1:
                    if len(obs) == 0:
                        obs.append(0)
                    if len(obs) > 1:
                        while len(obs) > 1:
                            obs.pop()
                for kss, vss in sondic1.items():
                    if "L2" in kss and sondic1[kss] != "         0.000  ":
                        obs.append(float(sondic1[kss][0:14]))
                if len(obs) != 2:
                    if len(obs) == 1:
                        obs.append(0)
                    if len(obs) > 2:
                        while len(obs) > 2:
                            obs.pop()
                for kss, vss in sondic1.items():
                    if "C1" in kss and sondic1[kss] != "         0.000  ":
                        obs.append(float(sondic1[kss][0:14]))
                    if "P1" in kss and sondic1[kss] != "         0.000  ":
                        obs.append(float(sondic1[kss][0:14]))
                if len(obs) != 3:
                    if len(obs) == 2:
                        obs.append(0)
                    if len(obs) > 3:
                        while len(obs) > 3:
                            obs.pop()
                for kss, vss in sondic1.items():
                    if "C2" in kss and sondic1[kss] != "         0.000  ":
                        obs.append(float(sondic1[kss][0:14]))
                    if "P2" in kss and sondic1[kss] != "         0.000  ":
                        obs.append(float(sondic1[kss][0:14]))
                if len(obs) != 4:
                    if len(obs) == 3:
                        obs.append(0)
                    if len(obs) > 4:
                        while len(obs) > 4:
                            obs.pop()
                lc[time_index.index(k)][prns.index(ks)][0] = obs[0]
                lc[time_index.index(k)][prns.index(ks)][1] = obs[1]
                lc[time_index.index(k)][prns.index(ks)][2] = obs[2]
                lc[time_index.index(k)][prns.index(ks)][3] = obs[3]
            if "S" in ks:
                sondic1 = sondic[ks]
                obs = []
                for kss, vss in sondic1.items():
                    if "L1" in kss and sondic1[kss] != "         0.000  ":
                        obs.append(float(sondic1[kss][0:14]))
                if len(obs) != 1:
                    if len(obs) == 0:
                        obs.append(0)
                    if len(obs) > 1:
                        while len(obs) > 1:
                            obs.pop()
                for kss, vss in sondic1.items():
                    if "L5" in kss and sondic1[kss] != "         0.000  ":
                        obs.append(float(sondic1[kss][0:14]))
                if len(obs) != 2:
                    if len(obs) == 1:
                        obs.append(0)
                    if len(obs) > 2:
                        while len(obs) > 2:
                            obs.pop()
                for kss, vss in sondic1.items():
                    if "C1" in kss and sondic1[kss] != "         0.000  ":
                        obs.append(float(sondic1[kss][0:14]))
                if len(obs) != 3:
                    if len(obs) == 2:
                        obs.append(0)
                    if len(obs) > 3:
                        while len(obs) > 3:
                            obs.pop()
                for kss, vss in sondic1.items():
                    if "C5" in kss and sondic1[kss] != "         0.000  ":
                        obs.append(float(sondic1[kss][0:14]))
                if len(obs) != 4:
                    if len(obs) == 3:
                        obs.append(0)
                    if len(obs) > 4:
                        while len(obs) > 4:
                            obs.pop()
                lc[time_index.index(k)][prns.index(ks)][0] = obs[0]
                lc[time_index.index(k)][prns.index(ks)][1] = obs[1]
                lc[time_index.index(k)][prns.index(ks)][2] = obs[2]
                lc[time_index.index(k)][prns.index(ks)][3] = obs[3]
            if "E" in ks:
                sondic1 = sondic[ks]
                obs = []
                for kss, vss in sondic1.items():
                    if "L1" in kss and sondic1[kss] != "         0.000  ":
                        obs.append(float(sondic1[kss][0:14]))
                if len(obs) != 1:
                    if len(obs) == 0:
                        obs.append(0)
                    if len(obs) > 1:
                        while len(obs) > 1:
                            obs.pop()
                for kss, vss in sondic1.items():
                    if "L5" in kss and sondic1[kss] != "         0.000  ":
                        obs.append(float(sondic1[kss][0:14]))
                if len(obs) != 2:
                    if len(obs) == 1:
                        obs.append(0)
                    if len(obs) > 2:
                        while len(obs) > 2:
                            obs.pop()
                for kss, vss in sondic1.items():
                    if "C1" in kss and sondic1[kss] != "         0.000  ":
                        obs.append(float(sondic1[kss][0:14]))
                if len(obs) != 3:
                    if len(obs) == 2:
                        obs.append(0)
                    if len(obs) > 3:
                        while len(obs) > 3:
                            obs.pop()
                for kss, vss in sondic1.items():
                    if "C5" in kss and sondic1[kss] != "         0.000  ":
                        obs.append(float(sondic1[kss][0:14]))
                if len(obs) != 4:
                    if len(obs) == 3:
                        obs.append(0)
                    if len(obs) > 4:
                        while len(obs) > 4:
                            obs.pop()
                lc[time_index.index(k)][prns.index(ks)][0] = obs[0]
                lc[time_index.index(k)][prns.index(ks)][1] = obs[1]
                lc[time_index.index(k)][prns.index(ks)][2] = obs[2]
                lc[time_index.index(k)][prns.index(ks)][3] = obs[3]
    return lc
def frelist(opath):
    prns = json.loads(oheader(opath))["PRNS"]
    Gf1 = 1.57542e9
    Gf2 = 1.2276e9
    c = 0.299792458e9
    Gl1 = c / Gf1
    Gl2 = c / Gf2
    Ga = (Gf1 / Gf2) ** 2
    Rslot = {"R01":	1,
    "R02":-4,
    "R03":5,
    "R04":6,
    "R05":1,
    "R06":-4,
    "R07":5,
    "R08":6,
    "R09":-2,
    "R10":-7,
    "R11":0,
    "R12":-1,
    "R13":-2,
    "R14":-7,
    "R15":0,
    "R16":-1,
    "R17":4,
    "R18":-3,
    "R19":3,
    "R20":2,
    "R21":4,
    "R22":-3,
    "R23":3,
    "R24":2
    }
    Ef1 = 1.57542e9
    Ef5 = 1.17645e9
    El1 = c / Ef1
    El5 = c / Ef5
    Ea = (Ef1 / Ef5) ** 2
    Sf1 = 1.57542e9
    Sf5 = 1.17645e9
    Sl1 = c / Sf1
    Sl5 = c / Sf5
    Sa = (Sf1 / Sf5) ** 2
    fres = np.zeros(shape=(len(prns), 5), dtype=float)
    for n in range(len(prns)):
        if prns[n][0] == "G":

            fres[n][0] = Gl1
            fres[n][1] = Gl2
            fres[n][2] = Ga
            fres[n][3] = Gf1**2/(Gf1**2-Gf2**2)
        if prns[n][0] == "R":
            for k,v in Rslot.items():
                if prns[n] == k:
                    Rf1 = (1602 + Rslot[k]*9/16)*10**6
                    Rf2 = (1246 + Rslot[k]*7/16)*10**6
                    fres[n][0] = c/Rf1
                    fres[n][1] = c/Rf2
                    fres[n][2] = (Rf1 / Rf2) ** 2
                    fres[n][3] = Rf1 ** 2 / (Rf1 ** 2 - Rf2 ** 2)
        if prns[n][0] == "E":

            fres[n][0] = El1
            fres[n][1] = El5
            fres[n][2] = Ea
            fres[n][3] = Ef1 ** 2 / (Ef1 ** 2 - Ef5 ** 2)
        if prns[n][0] == "S":

            fres[n][0] = Sl1
            fres[n][1] = Sl5
            fres[n][2] = Sa
            fres[n][3] = Sf1 ** 2 / (Sf1 ** 2 - Sf5 ** 2)
    return fres
def ION_MP(opath):
    odata = json.loads(observations(opath))
    prns = json.loads(oheader(opath))["PRNS"]
    maxoff = 10
    fres = frelist(opath)
    lc = LandC(opath)

    amount = len(prns)
    time_index = time_serises(odata)
    epoch_index = odata
    size = len(epoch_index)
    MP1MP2 = np.zeros(shape=(size, amount, 4), dtype=float)
    MP1mean = np.zeros(shape=(size, amount, 1), dtype=float)
    MP2mean = np.zeros(shape=(size, amount, 1), dtype=float)
    ini_svs = np.zeros(amount)
    for n in range(len(lc)):
        for m in range(len(lc[n])):
            if 0.0 in lc[n][m]:
                MP1MP2[n, m, 0] = 0
                MP1MP2[n, m, 1] = 0
            else:
                l1 = fres[m, 0]
                l2 = fres[m, 1]
                a = fres[m, 2]
                L1m = lc[n, m, 0] * l1
                L2m = lc[n, m, 1] * l2
                MP1MP2[n, m, 0] = lc[n, m, 2] - L1m*((2/(a-1))+1)+L2m*(2/(a-1))
                MP1MP2[n, m, 1] = lc[n, m, 3] - L1m*(2*a / (a - 1)) + L2m * (((2 * a) / (a - 1)) - 1)
                if ini_svs[m] == 0:
                    MP1MP2[n, m, 2] = (a / (a - 1)) * (L1m - L2m)
                    ini_svs[m] = MP1MP2[n, m, 2]
                    MP1MP2[n, m, 2] = 0
                else:
                    MP1MP2[n, m, 2] = (a / (a - 1)) * (L1m - L2m)
                    MP1MP2[n, m, 2] = ini_svs[m] - MP1MP2[n, m, 2]
                    if abs(MP1MP2[n, m, 2] > 10 * maxoff):
                        ini_svs[m] = (a / (a - 1)) * (L1m - L2m)
                    MP1MP2[n, m, 3] = ((MP1MP2[n, m, 2] - MP1MP2[n - 1, m, 2])
                                                            / (utf2gps(time_index[n].split())[1] -
                                                               utf2gps(time_index[n - 1].split())[
                                                                   1])) * 60
    for j in range(len(prns)):
        i = 0
        while i <= size and i!= size-1:
            while i <= size and i != size - 1:
                if MP1MP2[i][j][0] == 0:
                    i = i + 1
                elif abs(MP1MP2[i][j][0] - MP1MP2[i + 1][j][0]) < 2*maxoff:
                    l = i
                    ave1 = 0
                    ave2 = 0
                    while abs(MP1MP2[i][j][0] - MP1MP2[i + 1][j][0]) < 2*maxoff:
                        ave1 += MP1MP2[i][j][0]
                        ave2 += MP1MP2[i][j][1]
                        i = i + 1
                        if i >= size - 1:
                            break

                    ave1 = ave1 / (i - l)
                    ave2 = ave2 / (i - l)
                    for k in range(l, i+1):
                        if k != size:
                            MP1mean[k, j, 0] = ave1
                            MP2mean[k, j, 0] = ave2
                else:
                    MP1mean[i, j, 0] = 0
                    MP2mean[i, j, 0] = 0
                    i = i + 1
    for i in range(size):
        for j in range(len(prns)):
            if MP1MP2[i][j][0] != 0:
                MP1MP2[i][j][0] = MP1MP2[i][j][0] - MP1mean[i][j][0]
                MP1MP2[i][j][1] = MP1MP2[i][j][1] - MP2mean[i][j][0]
    epochs = time_index
    fieldnames = prns
    fieldnames.insert(0, "epoch")
    row1 = []
    row2 = []
    for n in range(len(MP1MP2)):
        element1 = {}
        element2 = {}
        element1["epoch"] = epochs[n]
        element2["epoch"] = epochs[n]
        for m in range(len(MP1MP2[n])):
            if MP1MP2[n][m][0] != 0.0 and MP1MP2[n][m][1] != 0.0:
                element1[prns[m + 1]] = [MP1MP2[n][m][0], MP1MP2[n][m][1]]
            if MP1MP2[n][m][2] != 0.0 and MP1MP2[n][m][3] != 0.0:
                element2[prns[m + 1]] = [MP1MP2[n][m][2], MP1MP2[n][m][3]]
        row1.append(element1)
        row2.append(element2)
    with open(opath[:-4] + "mp1mp2.csv", "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(row1)
    with open(opath[:-4] + "IonIod.csv", "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(row2)
    gps_prn = prns[1:]
    name = os.path.basename(opath)
    cloumn = math.ceil(len(prns)/18)
    plot(opath[:-4] + "MP1_plot.png", name[:-4] + " MP1 plot", gps_prn, epochs, MP1MP2, 0, cloumn, y_label="MP1 (m)")
    plot(opath[:-4] + "MP2_plot.png", name[:-4] + " MP2 plot", gps_prn, epochs, MP1MP2, 1, cloumn, y_label="MP2 (m)")
    plot(opath[:-4] + "ION_plot.png", name[:-4] + " ION plot", gps_prn, epochs, MP1MP2, 2, cloumn, y_label="ION (m)")
    plot(opath[:-4] + "IOD_plot.png", name[:-4] + " IOD plot", gps_prn, epochs, MP1MP2, 3, cloumn, y_label="IOD (m)")
    return MP1MP2
def cycleslip(opath):
    odata = json.loads(observations(opath))
    name = os.path.basename(opath)
    obs = LandC(opath)
    prns = json.loads(oheader(opath))["PRNS"]
    amount = len(prns)
    fres = frelist(opath)
    size = len(obs)
    PP = np.zeros(amount)
    PPhi = np.zeros(amount)
    P = np.zeros(shape=(size, amount), dtype=float)
    Phi = np.zeros(shape=(size, amount), dtype=float)
    deltaP = np.zeros(shape=(size, amount), dtype=float)
    deltaPhi = np.zeros(shape=(size, amount), dtype=float)
    cycleslips = np.zeros(shape=(size - 2, amount, 4), dtype=float)
    DP = np.zeros(shape=(size - 2, amount), dtype=float)
    DPhi = np.zeros(shape=(size - 2, amount), dtype=float)
    DDP = np.zeros(shape=(size - 4, amount), dtype=float)
    DDPhi = np.zeros(shape=(size - 4, amount), dtype=float)
    for n in range(size):
        for m in range(amount):
            if obs[n][m][0] !=0:
                if PP[m] == 0:
                    alpha1 = fres[m][3]
                    alpha2 = 1 - alpha1
                    lambda1 = fres[m][0]
                    lambda2 = fres[m][1]
                    P[n][m] = alpha1*obs[n][m][2] + alpha2*obs[n][m][3]
                    Phi[n][m] = alpha1*lambda1*obs[n][m][0]+alpha2*lambda2*obs[n][m][1]
                    PP[m] = P[n][m]
                    PPhi[m] = Phi[n][m]
                else:
                    P[n][m] = alpha1*obs[n][m][2] + alpha2*obs[n][m][3]
                    deltaP[n][m] = P[n][m] - PP[m]
                    Phi[n][m] = alpha1 * lambda1 * obs[n][m][0] + alpha2 * lambda2 * obs[n][m][1]
                    deltaPhi[n][m] = Phi[n][m] - PPhi[m]
    for n in range(1, size-1):
        for m in range(amount):
            cycleslips[n - 1][m][0] = (deltaP[n - 1][m] + deltaP[n + 1][m]) - 2 * deltaP[n][m]
            if abs(cycleslips[n-1][m][0]) > 280000:
                if cycleslips[n-1][m][0] > 0:
                    cycleslips[n-1][m][0] = cycleslips[n-1][m][0] - 299792.458
                else:
                    cycleslips[n - 1][m][0] = cycleslips[n - 1][m][0] + 299792.458
            cycleslips[n - 1][m][1] =(deltaPhi[n - 1][m] + deltaPhi[n + 1][m]) - 2 * deltaPhi[n][m]
    for n in range(1, size-3):
        for m in range(amount):
            cycleslips[n - 1][m][2] = (cycleslips[n - 1][m][0] + cycleslips[n + 1][m][0]) - 2 * cycleslips[n][m][0]
            cycleslips[n - 1][m][3] = (cycleslips[n - 1][m][1] + cycleslips[n + 1][m][1]) - 2 * cycleslips[n][m][1]

    epochs = time_serises(odata)
    fieldnames = prns
    fieldnames.insert(0, "epoch")
    row1 = []
    for n in range(len(cycleslips)):
        element1 = {}
        element1["epoch"] = epochs[n]
        for m in range(len(cycleslips[n])):
            if cycleslips[n][m][2] != 0.0 and cycleslips[n][m][3] != 0.0:
                element1[prns[m+1]] = cycleslips[n][m][3]
        row1.append(element1)
    with open(opath[:-4]+"CScarrier.csv", "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(row1)
    gps_prn = prns[1:]
    cloumn = math.ceil(len(prns) / 18)
    plot(opath[:-4] + "CycleSlipCarrier.png", name[:-4] + " Cycle Slip of Carrier-Phase observations plot", gps_prn,
         epochs, cycleslips, 3, cloumn,y_label="cycle slip (m)")
    return cycleslips
def QualityCheck(opath):
    # 获取文件后缀名
    _, ext = os.path.splitext(opath)

    # 构造npath
    npath = opath[:-1] + ('n' if ext[-1].lower() == 'o' else 'N') + ext[1:]

    # 检查npath是否存在
    if os.path.exists(npath):
        cyc = cycleslip(opath)
        ionmp = ION_MP(opath)
        aziele = azi_ele(opath)
        SatelliteSignalPlot(opath)
    else:
        cyc = cycleslip(opath)
        ionmp = ION_MP(opath)
        SatelliteSignalPlot(opath)
def batchQC(rootpath, keywords_list, extension):
    PINEXS = DataFinding(rootpath, keywords_list, extension)
    for RINEX in RINEXS:
        QualityCheck(RINEX)


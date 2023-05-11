import os
import re
import datetime
from PyRINEX.bookworm import detectCode as dc
from PyRINEX.bookworm import readtxt as rt
from PyRINEX.bookworm import file_name_extesion_judgement as fnej
from PyRINEX.bookworm import file_name_extesion_judgement_file as fnejf
from PyRINEX.bookworm import show_files as sf
from PyRINEX.bookworm import readtxt_in_encoding as rte
from PyRINEX.bookworm import readtxt as rt
from PyRINEX.venom import RinexO as RO
from PyRINEX.venom import rinex2_obs as r2o
from PyRINEX.venom import rinex3_obs as r3o
from PyRINEX.venom import rinex2_nav as r2n
from PyRINEX.venom import rinex3_nav as r3n
import numpy as np
import math
from numpy import *
import csv
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from dateutil import parser
import time

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
def SatelliteSignalPlot(path):
    name = os.path.basename(path)
    version = RO(path).version()
    obss = ""
    if version == "2":
        obss = r2o(path).epoch_index()
    if version == "3":
        obss = r3o(path).epoch_index()
    id = []
    TIME = []
    for m in range(len(obss)):
        for k, v in obss[m].items():
            if type(obss[m][k]) == dict:
                TIME.append(":".join(obss[m]["time"].split()[3:6]))
                id.append(k)
    ax = plt.subplot()
    hours = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(hours)
    ax.set_title(name + " Signal Plot")
    time = [parser.parse(x) for x in TIME]
    plt.scatter(time, id, s=2, marker=".", alpha=0.75, c="g")
    plt.grid(linestyle = '--', linewidth=0.5)
    plt.ylabel("SATELLITE NO")
    plt.gcf().set_size_inches(15, 7)
    filename = path[:-4] + "SignalPlot.jpg"
    plt.savefig(filename, dpi=360)

    plt.cla()
    plt.clf()
def time_serises(path):
    version = RO(path).version()
    if version == "2":
        obss = r2o(path).epoch_index()
    if version == "3":
        obss = r3o(path).epoch_index()
    list = []
    for n in range(len(obss)):
        list.append(obss[n]["time"])
    return list
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
def gps_pos(path):
    pos_list = []
    a = 6378137
    f = 1/298.257223563
    Omega_e_Dot = 7.2921151467e-5
    mu = 3.986005e14
    c = 2.99792458e8
    version = RO(path).version()
    if path[-1] == "O":

        nfile_path = path[0:-1] + "N"
    if path[-1] == "o":

        nfile_path = path[0:-1]+"n"
    if version == "2":
        epoch_index = r2o(path).epoch_index()
        Nlist = r2n(nfile_path).sat_nav()
    if version == "3":
        epoch_index = r3o(path).epoch_index()
        Nlist = r3n(nfile_path).sat_nav()
    postions = []
    for n in range(len(epoch_index)):
        postions_per_epoch = {}
        postions_per_epoch["epoch"] = epoch_index[n]["time"]
        gps = []
        for k, v in epoch_index[n].items():
            utf = epoch_index[n]["time"].split()
            gps_seconds = utf2gps(utf)[1]
            if "G" in k:
                gps.append(k)
        for m in range(len(gps)):
            index = []
            t_kResult = []
            for l in range(len(Nlist)):
                for k, v in Nlist[l].items():
                    gps_prn = k
                    dic = Nlist[l][gps_prn]
                toe = float(dic["Toe Time of Ephemeris"].replace("D", "E"))
                try:
                    if gps[m] == gps_prn and t_kCal(gps_seconds, toe) != "not valid":

                        index.append(l)
                        t_kResult.append(abs(t_kCal(gps_seconds, toe)))
                except ValueError:
                    if gps[m] == gps_prn:
                        index.append(l)
                        t_kResult.append(abs(t_kCal(gps_seconds, toe)))

            dic = Nlist[index[t_kResult.index(min(t_kResult))]][gps[m]]
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
            postions_per_epoch[gps[m]] = [x_k, y_k, z_k]
        postions.append(postions_per_epoch)
    return postions
def azi_ele(path):
    name = os.path.basename(path)
    filename = path[:-4] + "SkyPlot.jpg"
    xyz = RO(path).xyz()
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    p = math.sqrt(x ** 2 + y ** 2)
    R = math.sqrt(x ** 2 + y ** 2 + z ** 2)
    e = [[-y / p, x / p, 0],
        [-x * z / (p * R), -y * z / (p * R), p / R],
        [x / R, y / R, z / R]]
    postion = gps_pos(path)
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
    ax = plt.subplot(111, polar=True)
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N')
    ax.set_rlim(90, 0)
    ax.set_title(name+"SkyPlot")
    element = {}
    version=RO(path).version()
    if version == "2":
        gpss = r2o(path).gpslist()
    if version == "3":
        gpss = r3o(path).gpslist()
    colormap = plt.cm.nipy_spectral
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
    plt.savefig(filename, dpi=360, bbox_inches='tight')
    plt.cla()
    plt.clf()
    fieldnames = gpss
    print(fieldnames)
    fieldnames.insert(0, "epoch")
    row = []
    for n in range(len(azi_ele_list)):
        element = {}
        element["epoch"] = azi_ele_list[n]["epoch"]
        for k, v in azi_ele_list[n].items():
            element[k] = azi_ele_list[n][k]
        row.append(element)
    with open(path[:-4]+"aziele.csv", "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames = fieldnames)
        writer.writeheader()
        writer.writerows(row)
    return azi_ele_list
def plot(filename, title, gpsserises, epochs, dataset, type):
    TIME = []
    for m in range(len(epochs)):
        TIME.append(":".join(epochs[m].split()[3:6]))
    ax = plt.subplot()
    hours = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(hours)
    ax.set_title(title)
    time = [parser.parse(x) for x in TIME]
    element = {}
    gps_prn = gpsserises
    number_of_plots = len(gps_prn)
    colormap = plt.cm.nipy_spectral
    colors = [colormap(i) for i in np.linspace(0, 1, number_of_plots)]
    ax.set_prop_cycle('color', colors)
    for n in range(len(gps_prn)):
        element[gps_prn[n]] = []
    for n in range(len(dataset)):
        for m in range(len(gps_prn)):
            if dataset[n][m][type] != 0 and abs(dataset[n][m][type]) < 100:
            #if dataset[n][m][type] != 0:
                element[gps_prn[m]].append([time[n], dataset[n][m][type]])
    for n in range(len(gps_prn)):
        data = element[gps_prn[n]]
        x = []
        y = []
        for m in range(len(data)):
            x.append(data[m][0])
            y.append(data[m][1])
        ax.scatter(x, y, s=1, label=gps_prn[n])
    plt.legend(bbox_to_anchor=(1.05, 0), loc=3, borderaxespad=0)
    plt.savefig(filename, dpi=360, bbox_inches='tight')
    plt.cla()
def ion(path):
    aziele = azi_ele(path)
    version = RO(path).version()
    c = 2.99792458e8
    if path[-1] == "O":

        nfile_path = path[0:-1] + "N"
    if path[-1] == "o":

        nfile_path = path[0:-1]+"n"
    if version == "2":
        epoch_index = r2o(path).epoch_index()
        ionpara = r2n(nfile_path).IONpar()
    if version == "3":
        epoch_index = r3o(path).epoch_index()
        ionpara = r3n(nfile_path).IONpar()
    if ionpara[0] == "NO ION data.":
        ionpara = ["0.1118D-07", "-0.7451D-08", "-0.5961D-07", "0.1192D-06",
	    "0.1167D+06", "-0.2294D+06", "-0.1311E+06", "0.1049E+07"]
    a = ionpara[0:4]
    b = ionpara[4:8]
    alfa = []
    beta = []
    for n in range(len(a)):
        alfa.append(float(a[n].replace("D", "E")))
        beta.append(float(b[n].replace("D", "E")))
    coordinate = xyz2BLH(RO(path).xyz()[0], RO(path).xyz()[1], RO(path).xyz()[2],
                         rad=True)
    lon, lat = coordinate[0], coordinate[1]
    for n in range(len(aziele)):
        epoch = aziele[n]["epoch"].split()
        td = utf2gps(epoch)[0]
        gps = []
        for k, v in aziele[n].items():
            if "G" in k:
                gps.append(k)
        for m in range(len(gps)):
            azi = (aziele[n][gps[m]][0])
            ele = (aziele[n][gps[m]][1])
            psi = 0.0137/(ele/math.pi+0.11) - 0.022
            phi = lat/math.pi+psi*math.cos(azi)
            if phi > 0.416:
                phi = 0.416
            if phi < -0.416:
                phi = -0.416
            lam = lon/math.pi + psi*math.sin(azi)/math.cos(phi*math.pi)
            phi=phi+0.064*math.cos((lam-1.617)*math.pi)
            tt = lam * 43200 + td
            tt = tt-math.floor(tt/86400)*86400
            f = 1 + 16*((0.53-ele/math.pi)**3)
            amp = alfa[0]+phi*(alfa[1]+phi*(alfa[2]+phi*alfa[3]))
            per = beta[0]+phi*(beta[1]+phi*(beta[2]+phi*beta[3]))
            if amp < 0:
                amp = 0
            if per < 72000:
                per = 72000
            x = 2*math.pi*(tt-50400)/per
            if abs(x) < 1.57:
                ion = c*f*(5*10**(-9)+amp*(1-(x**2/2)+(x**4/24)))
            else:
                ion = c*f(5*10**(-9))
            print(gps[m], ion)
def DualIon_mp(path):
    name = os.path.basename(path)
    f1 = 1.57542e9
    f2 = 1.2276e9
    c = 0.299792458e9
    l1 = c / f1
    l2 = c / f2
    a = (f1 / f2)**2
    maxoff = 10
    version = RO(path).version()
    allresult = []
    if version == "2":
        epoch_dictionary = r2o(path).epoch_index()
        size = len(epoch_dictionary)
        gpss = r2o(path).gpslist()
        gpsamount = len(gpss)
        MP1MP2 = np.zeros(shape=(size, gpsamount, 4), dtype=float)
        MP1mean = np.zeros(shape=(size, gpsamount, 1), dtype=float)
        MP2mean = np.zeros(shape=(size, gpsamount, 1), dtype=float)
        ini_svs = np.zeros(gpsamount)
        for n in range(len(epoch_dictionary)):
            gps_prn = []
            for k, v in epoch_dictionary[n].items():
                epoch = epoch_dictionary[n]["time"].split()
                if "G" in k:
                    gps_prn.append(k)
            for m in range(len(gps_prn)):
                row_number = int(gps_prn[m][1:3])
                obstypes = epoch_dictionary[n][gps_prn[m]]
                too = [] #types of observ
                for k, v in obstypes.items():
                    too.append(k)
                obs = []
                obs.append(float(obstypes["L1"][0:14]))
                obs.append(float(obstypes["L2"][0:14]))
                if "P1" in too:
                    if obstypes["P1"][9:14] == "0.000":
                        if "C1" in too:
                            obs.append(float(obstypes["C1"][0:14]))
                        else:
                            obs.append(float(obstypes["P1"][0:14]))
                    else:
                        obs.append(float(obstypes["P1"][0:14]))
                else:
                    obs.append(float(obstypes["C1"][0:14]))
                if "P2" in too:
                    if obstypes["P2"][9:14] == "0.000":
                        if "C2" in too:
                            obs.append(float(obstypes["C2"][0:14]))
                        else:
                            obs.append(float(obstypes["P2"][0:14]))
                    else:
                        obs.append(float(obstypes["P2"][0:14]))
                else:
                    obs.append(float(obstypes["C2"][0:14]))
                if 0.0 in obs:
                    MP1MP2[n, gpss.index(gps_prn[m]), 0] = 0
                    MP1MP2[n, gpss.index(gps_prn[m]), 1] = 0
                else:
                    L1m = obs[0]*l1
                    L2m = obs[1]*l2
                    MP1MP2[n, gpss.index(gps_prn[m]), 0] = obs[2]-L1m*((2/(a-1))+1)+L2m*(2/(a-1))
                    MP1MP2[n, gpss.index(gps_prn[m]), 1] = obs[3]-L1m*(2*a/(a-1))+L2m*(((2*a)/(a-1))-1)
                    if ini_svs[gpss.index(gps_prn[m])]==0:
                        MP1MP2[n, gpss.index(gps_prn[m]), 2] = (a/(a-1))*(L1m-L2m)
                        ini_svs[gpss.index(gps_prn[m])] = MP1MP2[n, gpss.index(gps_prn[m]), 2]
                        MP1MP2[n, gpss.index(gps_prn[m]), 2] = 0
                    else:
                        MP1MP2[n, gpss.index(gps_prn[m]), 2] = (a / (a - 1)) * (L1m - L2m)
                        MP1MP2[n, gpss.index(gps_prn[m]), 2] = ini_svs[gpss.index(gps_prn[m])]-MP1MP2[n, gpss.index(gps_prn[m]), 2]
                        if abs(MP1MP2[n, gpss.index(gps_prn[m]), 2] > 10*maxoff):
                            ini_svs[gpss.index(gps_prn[m])] = (a / (a - 1)) * (L1m - L2m)
                        MP1MP2[n, gpss.index(gps_prn[m]), 3] = ((MP1MP2[n, gpss.index(gps_prn[m]), 2]-
                                                                MP1MP2[n-1, gpss.index(gps_prn[m]), 2])
                        /(utf2gps(epoch_dictionary[n]["time"].split())[1]-utf2gps(epoch_dictionary[n-1]["time"].split())[1]))*60

        for j in range(len(gpss)):
            i = 0
            while i <= size and i!= size-1:
                if MP1MP2[i][j][0] == 0:
                    i = i+1
                elif abs(MP1MP2[i][j][0] - MP1MP2[i+1][j][0]) < 10:
                    l = i
                    ave1 = 0
                    ave2 = 0
                    ave1 = ave1 + MP1MP2[i][j][0]
                    ave2 = ave2 + MP1MP2[i][j][1]
                    i = i+1
                    if i < size - 1:
                        while abs(MP1MP2[i][j][0] - MP1MP2[i+1][j][0]) < 10:
                            ave1 = ave1 + MP1MP2[i][j][0]
                            ave2 = ave2 + MP1MP2[i][j][1]
                            i = i+1
                            if i >= size -1:
                                break
                    #while abs(MP1MP2[i][j][0] - MP1MP2[i + 1][j][0]) < 10:
                    #    ave1 = ave1 + MP1MP2[i][j][0]
                    #    ave2 = ave2 + MP1MP2[i][j][1]
                    #    i = i+1
                    #    if i >= size -1:
                    #        break
                    ave1 = ave1/(i-l)
                    ave2 = ave2 / (i - l)
                    for k in range(l,i):
                        #if k != size:
                        MP1mean[k, j, 0] = ave1
                        #print(MP1mean[i, j, 0])
                        MP2mean[k, j, 0] = ave2

                else:
                    MP1mean[i,j,0] = 0
                    MP2mean[i,j,0] = 0
                    i =i+1
        for i in range(size):
            for j in range(len(gpss)):
                if MP1MP2[i][j][0] != 0:

                    MP1MP2[i][j][0] = MP1MP2[i][j][0] - MP1mean[i][j][0]
                    MP1MP2[i][j][1] = MP1MP2[i][j][1] - MP2mean[i][j][0]
    if version == "3":
        epoch_dictionary = r3o(path).epoch_index()
        size = len(epoch_dictionary)
        gpss = r3o(path).gpslist()
        gpsamount = len(gpss)
        MP1MP2 = np.zeros(shape=(size, gpsamount, 4), dtype=float)
        MP1mean = np.zeros(shape=(size, gpsamount, 1), dtype=float)
        MP2mean = np.zeros(shape=(size, gpsamount, 1), dtype=float)
        ini_svs = np.zeros(gpsamount)
        for n in range(len(epoch_dictionary)):
            gps_prn = []
            for k, v in epoch_dictionary[n].items():
                epoch = epoch_dictionary[n]["time"]
                if "G" in k:
                    gps_prn.append(k)
            for m in range(len(gps_prn)):
                row_number = int(gps_prn[m][1:3])
                obstypes = epoch_dictionary[n][gps_prn[m]]
                obs = []
                for k, v in obstypes.items():
                    if "L1" in k and obstypes[k] != "         0.000  ":
                       obs.append(float(obstypes[k][0:14]))
                if len(obs)!=1:
                    if len(obs) == 0:
                        obs.append(0)
                    if len(obs) > 1:
                        while len(obs) > 1:
                            obs.pop()
                for k, v in obstypes.items():
                    if "L2" in k and obstypes[k] != "         0.000  ":
                       obs.append(float(obstypes[k][0:14]))
                if len(obs)!=2:
                    if len(obs) == 1:
                        obs.append(0)
                    if len(obs) > 2:
                        while len(obs) > 2:
                            obs.pop()
                for k, v in obstypes.items():
                    if "C1" in k and obstypes[k] != "         0.000  ":
                       obs.append(float(obstypes[k][0:14]))
                if len(obs)!=3:
                    if len(obs) == 2:
                        obs.append(0)
                    if len(obs) > 3:
                        while len(obs) > 3:
                            obs.pop()
                for k, v in obstypes.items():
                    if "C2" in k and obstypes[k] != "         0.000  ":
                       obs.append(float(obstypes[k][0:14]))
                if len(obs)!=4:
                    if len(obs) == 3:
                        obs.append(0)
                    if len(obs) > 4:
                        while len(obs) > 4:
                            obs.pop()
                if 0.0 in obs:
                    MP1MP2[n, gpss.index(gps_prn[m]), 0] = 0
                    MP1MP2[n, gpss.index(gps_prn[m]), 1] = 0
                else:
                    L1m = obs[0]*l1
                    L2m = obs[1]*l2
                    MP1MP2[n, gpss.index(gps_prn[m]), 0] = obs[2]-L1m*((2/(a-1))+1)+L2m*(2/(a-1))
                    MP1MP2[n, gpss.index(gps_prn[m]), 1] = obs[3]-L1m*(2*a/(a-1))+L2m*(((2*a)/(a-1))-1)
                    if ini_svs[gpss.index(gps_prn[m])] == 0:
                        MP1MP2[n, gpss.index(gps_prn[m]), 2] = (a / (a - 1)) * (L1m - L2m)
                        ini_svs[gpss.index(gps_prn[m])] = MP1MP2[n, gpss.index(gps_prn[m]), 2]
                        MP1MP2[n, gpss.index(gps_prn[m]), 2] = 0
                    else:
                        MP1MP2[n, gpss.index(gps_prn[m]), 2] = (a / (a - 1)) * (L1m - L2m)
                        MP1MP2[n, gpss.index(gps_prn[m]), 2] = ini_svs[gpss.index(gps_prn[m])] - MP1MP2[
                            n, gpss.index(gps_prn[m]), 2]
                        if abs(MP1MP2[n, gpss.index(gps_prn[m]), 2] > 10 * maxoff):
                            ini_svs[gpss.index(gps_prn[m])] = (a / (a - 1)) * (L1m - L2m)
                        MP1MP2[n, gpss.index(gps_prn[m]), 3] = ((MP1MP2[n, gpss.index(gps_prn[m]), 2] -
                                                                 MP1MP2[n - 1, gpss.index(gps_prn[m]), 2])
                                                                / (utf2gps(epoch_dictionary[n]["time"].split())[1] -
                                                                   utf2gps(epoch_dictionary[n - 1]["time"].split())[
                                                                       1])) * 60
        for j in range(len(gpss)):
            i = 0
            while i <= size and i!= size-1:
                if MP1MP2[i][j][0] == 0:
                    i = i+1
                elif abs(MP1MP2[i][j][0] - MP1MP2[i+1][j][0]) < 10:
                    l = i
                    ave1 = 0
                    ave2 = 0
                    ave1 = ave1 + MP1MP2[i][j][0]
                    ave2 = ave2 + MP1MP2[i][j][1]
                    i = i+1
                    if i < size-1:
                        while abs(MP1MP2[i][j][0] - MP1MP2[i+1][j][0]) < 10:
                            ave1 = ave1 + MP1MP2[i][j][0]
                            ave2 = ave2 + MP1MP2[i][j][1]
                            i = i+1

                            if i >= size -1:
                                break
                    ave1 = ave1/(i-l)
                    ave2 = ave2 / (i - l)
                    for k in range(l,i):
                        #if k != size:
                        MP1mean[k, j, 0] = ave1
                        MP2mean[k, j, 0] = ave2

                else:
                    MP1mean[i,j,0] = 0
                    MP2mean[i,j,0] = 0
                    i =i+1
        for i in range(size):
            for j in range(len(gpss)):
                if MP1MP2[i][j][0] != 0:
                    #print(MP1MP2[i][j][0] , MP1mean[i][j][0])
                    MP1MP2[i][j][0] = MP1MP2[i][j][0] -MP1mean[i][j][0]
                    MP1MP2[i][j][1] = MP1MP2[i][j][1] - MP2mean[i][j][0]
    MP1MP2 = MP1MP2[:-1]
    epochs = time_serises(path)
    if version == "2":
        gpss = r2o(path).gpslist()
    if version == "3":
        gpss = r3o(path).gpslist()
    fieldnames = gpss
    fieldnames.insert(0, "epoch")
    row1 = []
    row2 = []
    for n in range(len(MP1MP2)):
        element1 = {}
        element2 = {}
        element1["epoch"] = epochs[n]
        element2["epoch"] = epochs[n]
        for m in range(len(MP1MP2[n])):
            if MP1MP2[n][m][0] !=0.0 and MP1MP2[n][m][1] !=0.0:
                element1[gpss[m+1]] = [MP1MP2[n][m][0], MP1MP2[n][m][1]]
            if MP1MP2[n][m][2] != 0.0 and MP1MP2[n][m][3] != 0.0:
                element2[gpss[m+1]] = [MP1MP2[n][m][2], MP1MP2[n][m][3]]
        row1.append(element1)
        row2.append(element2)
    with open(path[:-4]+"mp1mp2.csv", "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(row1)
    with open(path[:-4]+"IonIod.csv", "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(row2)
    gps_prn = gpss[1:]
    plot(path[:-4] + "MP1_plot.jpg", name + " MP1 plot", gps_prn, epochs, MP1MP2, 0)
    plot(path[:-4] + "MP2_plot.jpg", name + " MP2 plot", gps_prn, epochs, MP1MP2, 1)
    plot(path[:-4] + "ION_plot.jpg", name + " ION plot", gps_prn, epochs, MP1MP2, 2)
    plot(path[:-4] + "IOD_plot.jpg", name + " IOD plot", gps_prn, epochs, MP1MP2, 3)
    return MP1MP2
def LAndC(path):
    version = RO(path).version()
    if version == "2":
        epoch_dictionary = r2o(path).epoch_index()
        size = len(epoch_dictionary)
        gpss = r2o(path).gpslist()
        gpsamount = len(gpss)
        lc = np.zeros(shape=(size, gpsamount, 4), dtype=float)
        for n in range(len(epoch_dictionary)):
            gps_prn = []
            for k, v in epoch_dictionary[n].items():
                epoch = epoch_dictionary[n]["time"]
                if "G" in k:
                    gps_prn.append(k)
            for m in range(len(gps_prn)):
                row_number = int(gps_prn[m][1:3])
                obstypes = epoch_dictionary[n][gps_prn[m]]
                too = []  # types of observ
                for k, v in obstypes.items():
                    too.append(k)
                obs = []
                obs.append(float(obstypes["L1"][0:14]))
                obs.append(float(obstypes["L2"][0:14]))
                if "P1" in too:
                    if obstypes["P1"][9:14] == "0.000":
                        if "C1" in too:
                            obs.append(float(obstypes["C1"][0:14]))
                        else:
                            obs.append(float(obstypes["P1"][0:14]))
                    else:
                        obs.append(float(obstypes["P1"][0:14]))
                else:
                    obs.append(float(obstypes["C1"][0:14]))
                if "P2" in too:
                    if obstypes["P2"][9:14] == "0.000":
                        if "C2" in too:
                            obs.append(float(obstypes["C2"][0:14]))
                        else:
                            obs.append(float(obstypes["P2"][0:14]))
                    else:
                        obs.append(float(obstypes["P2"][0:14]))
                else:
                    obs.append(float(obstypes["C2"][0:14]))
                lc[n][gpss.index(gps_prn[m])][0] = obs[0]
                lc[n][gpss.index(gps_prn[m])][1] = obs[1]
                lc[n][gpss.index(gps_prn[m])][2] = obs[2]
                lc[n][gpss.index(gps_prn[m])][3] = obs[3]
        return lc
    if version == "3":
        epoch_dictionary = r3o(path).epoch_index()
        size = len(epoch_dictionary)
        gpss = r3o(path).gpslist()
        gpsamount = len(gpss)
        lc = np.zeros(shape=(size, gpsamount, 4), dtype=float)
        for n in range(len(epoch_dictionary)):
            gps_prn = []
            for k, v in epoch_dictionary[n].items():
                epoch = epoch_dictionary[n]["time"]
                if "G" in k:
                    gps_prn.append(k)
            for m in range(len(gps_prn)):
                row_number = int(gps_prn[m][1:3])
                obstypes = epoch_dictionary[n][gps_prn[m]]
                obs = []
                for k, v in obstypes.items():
                    if "L1" in k and obstypes[k] != "         0.000  ":
                        obs.append(float(obstypes[k][0:14]))
                if len(obs) != 1:
                    if len(obs) == 0:
                        obs.append(0)
                    if len(obs) > 1:
                        while len(obs) > 1:
                            obs.pop()
                for k, v in obstypes.items():
                    if "L2" in k and obstypes[k] != "         0.000  ":
                        obs.append(float(obstypes[k][0:14]))
                if len(obs) != 2:
                    if len(obs) == 1:
                        obs.append(0)
                    if len(obs) > 2:
                        while len(obs) > 2:
                            obs.pop()
                for k, v in obstypes.items():
                    if "C1" in k and obstypes[k] != "         0.000  ":
                        obs.append(float(obstypes[k][0:14]))
                if len(obs) != 3:
                    if len(obs) == 2:
                        obs.append(0)
                    if len(obs) > 3:
                        while len(obs) > 3:
                            obs.pop()
                for k, v in obstypes.items():
                    if "C2" in k and obstypes[k] != "         0.000  ":
                        obs.append(float(obstypes[k][0:14]))
                if len(obs) != 4:
                    if len(obs) == 3:
                        obs.append(0)
                    if len(obs) > 4:
                        while len(obs) > 4:
                            obs.pop()
                lc[n][gpss.index(gps_prn[m])][0] = obs[0]
                lc[n][gpss.index(gps_prn[m])][1] = obs[1]
                lc[n][gpss.index(gps_prn[m])][2] = obs[2]
                lc[n][gpss.index(gps_prn[m])][3] = obs[3]
        return lc
def cycleslip(path):
    c = 0.299792458e9
    f0 = 10.23*1000000
    f1 = 154*f0
    f2 = 120*f0
    lambda1 = c/f1
    lambda2 = c/f2
    alpha1 = f1**2/(f1**2-f2**2)
    alpha2 = 1-alpha1
    obs = LAndC(path)
    gpsamount = len(obs[0])
    size = len(obs)
    PP = np.zeros(gpsamount)
    PPhi = np.zeros(gpsamount)
    P = np.zeros(shape=(size, gpsamount), dtype=float)
    Phi = np.zeros(shape=(size, gpsamount), dtype=float)
    deltaP = np.zeros(shape=(size, gpsamount), dtype=float)
    deltaPhi = np.zeros(shape=(size, gpsamount), dtype=float)
    cycleslips = np.zeros(shape=(size - 2, gpsamount, 4), dtype=float)
    DP = np.zeros(shape=(size - 2, gpsamount), dtype=float)
    DPhi = np.zeros(shape=(size - 2, gpsamount), dtype=float)
    DDP = np.zeros(shape=(size - 4, gpsamount), dtype=float)
    DDPhi = np.zeros(shape=(size - 4, gpsamount), dtype=float)
    for n in range(size):
        start = time.time()
        for m in range(gpsamount):
            if obs[n][m][0] !=0:
                if PP[m] == 0:
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
            for m in range(gpsamount):
                cycleslips[n - 1][m][0] = (deltaP[n - 1][m] + deltaP[n + 1][m]) - 2 * deltaP[n][m]
                if abs(cycleslips[n-1][m][0]) > 280000:
                    if cycleslips[n-1][m][0] > 0:
                        cycleslips[n-1][m][0] = cycleslips[n-1][m][0] - 299792.458
                    else:
                        cycleslips[n - 1][m][0] = cycleslips[n - 1][m][0] + 299792.458
                cycleslips[n - 1][m][1] =(deltaPhi[n - 1][m] + deltaPhi[n + 1][m]) - 2 * deltaPhi[n][m]
        for n in range(1, size-3):
            for m in range(gpsamount):
                cycleslips[n - 1][m][2] = (cycleslips[n - 1][m][0] + cycleslips[n + 1][m][0]) - 2 * cycleslips[n][m][0]
                cycleslips[n - 1][m][3] = (cycleslips[n - 1][m][1] + cycleslips[n + 1][m][1]) - 2 * cycleslips[n][m][1]
        end = time.time()
    epochs = time_serises(path)
    version = RO(path).version()
    if version == "2":
        gpss = r2o(path).gpslist()
    if version == "3":
        gpss = r3o(path).gpslist()
    fieldnames = gpss
    fieldnames.insert(0, "epoch")
    row1 = []
    for n in range(len(cycleslips)):
        element1 = {}
        element1["epoch"] = epochs[n]
        for m in range(len(cycleslips[n])):
            if cycleslips[n][m][2] != 0.0 and cycleslips[n][m][3] != 0.0:
                #element1[gpss[m+1]] = [cycleslips[n][m][2], cycleslips[n][m][3]]
                element1[gpss[m + 1]] = cycleslips[n][m][3]
        row1.append(element1)
    with open(path+"cyc.csv", "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(row1)
    gps_prn = gpss[1:]
    name = os.path.basename(path)
    plot(path + "cyc.jpg", name + " Cycle Slip of Carrier-Phase observations plot", gps_prn, epochs, cycleslips, 3)
    return cycleslips
def QualityCheck(paths):
    for n in range(len(paths)):
        path = paths[n]
        name = os.path.basename(paths[n])
        version = RO(paths[n]).version()
        if version == "2":
            gpss = r2o(paths[n]).gpslist()
        if version == "3":
            gpss = r3o(paths[n]).gpslist()
        ele_dict = {}
        MP1_dict = {}
        MP2_dict = {}
        ION_dict = {}
        IOD_dict = {}
        cyc_dict = {}
        for m in range(len(gpss)):
            ele_dict[gpss[m]] = []
            MP1_dict[gpss[m]] = []
            MP2_dict[gpss[m]] = []
            ION_dict[gpss[m]] = []
            IOD_dict[gpss[m]] = []
            cyc_dict[gpss[m]] = []
        SatelliteSignalPlot(paths[n])
        aziele = azi_ele(paths[n])
        for m in range(len(aziele)):
            for k, v in aziele[m].items():
                for keys, values in ele_dict.items():
                    if k == keys:
                        ele_dict[k].append(aziele[m][k][1])
        MP1MP2 = DualIon_mp(paths[n])
        cyc = cycleslip(paths[n])
        for m in range(len(MP1MP2)):
            for o in range(len(gpss)):
                if MP1MP2[m][o][0] != 0:
                    MP1_dict[gpss[o]].append(MP1MP2[m][o][0])
                if MP1MP2[m][o][1] != 0:
                    MP2_dict[gpss[o]].append(MP1MP2[m][o][1])
                if MP1MP2[m][o][2] != 0:
                    ION_dict[gpss[o]].append(MP1MP2[m][o][2])
                if MP1MP2[m][o][3] != 0:
                    IOD_dict[gpss[o]].append(MP1MP2[m][o][3])
        for m in range(len(cyc)):
            for o in range(len(gpss)):
                if cyc[m][o][3] != 0:
                    cyc_dict[gpss[o]].append(cyc[m][o][3])
        ele_check = {}
        ele_check["TYPE"] = "elevation result"
        allvalues = []
        allunpassed = []
        for k, v in ele_dict.items():
            for m in range(len(gpss)):
                amount = len(ele_dict[gpss[m]])
                unpassed = []

                for n in range(len(ele_dict[gpss[m]])):
                    if ele_dict[gpss[m]][n] < 0.2617994:
                        unpassed.append(ele_dict[gpss[m]][n])
                        allunpassed.append(ele_dict[gpss[m]][n])
                    allvalues.append(ele_dict[gpss[m]][n])
                ele_check[gpss[m]] = 1 - len(unpassed)/amount
        ele_check["TOTAL"] = 1-len(allunpassed)/len(allvalues)
        #for k, v in ele_dict.items():
        #    for m in range(len(gpss)):
        #        for n in range(len(ele_dict[gpss[m]])):
        #            ele_dict[gpss[m]][n] = abs(ele_dict[gpss[m]][n])
        #        ele_check[gpss[m]] = mean(ele_dict[gpss[m]][n])
        MP1_check = {}
        MP1_check["TYPE"] = "MP1 result"
        MP1_values1 = []
        MP1_values2 = []
        for m in range(len(gpss)):
            amount = len(MP1_dict[gpss[m]])
            unpassed = []
            for n in range(len(MP1_dict[gpss[m]])):
                if abs(MP1_dict[gpss[m]][n]) > 1:
                    unpassed.append(MP1_dict[gpss[m]][n])
                    MP1_values1.append(MP1_dict[gpss[m]][n])
                MP1_values2.append(MP1_dict[gpss[m]][n])
            MP1_check[gpss[m]] = 1 - len(unpassed)/amount
        MP1_check["TOTAL"] = 1-len(MP1_values1)/len(MP1_values2)
        #MP1_values = []
#
#
        #for m in range(len(gpss)):
        #    for n in range(len(MP1_dict[gpss[m]])):
        #        MP1_dict[gpss[m]][n] = abs(MP1_dict[gpss[m]][n])
        #        MP1_values.append(abs(MP1_dict[gpss[m]][n]))
        #    MP1_check[gpss[m]] = mean(MP1_dict[gpss[m]])
        #MP1_check["TOTAL"] = mean(MP1_values)
        MP2_check = {}
        MP2_check["TYPE"] = "MP2 result"
        MP2_values1 = []
        MP2_values2 = []
        for m in range(len(gpss)):
            amount = len(MP2_dict[gpss[m]])
            unpassed = []
            for n in range(len(MP2_dict[gpss[m]])):
                if abs(MP2_dict[gpss[m]][n]) > 2:
                    unpassed.append(MP2_dict[gpss[m]][n])
                    MP2_values1.append(MP2_dict[gpss[m]][n])
                MP2_values2.append(MP2_dict[gpss[m]][n])
            MP2_check[gpss[m]] = 1 - len(unpassed)/amount
        MP2_check["TOTAL"] = 1-len(MP2_values1)/len(MP2_values2)
        #MP2_values = []
#
#
        #for m in range(len(gpss)):
        #    for n in range(len(MP1_dict[gpss[m]])):
        #        MP2_dict[gpss[m]][n] = abs(MP2_dict[gpss[m]][n])
        #        MP2_values.append(abs(MP2_dict[gpss[m]][n]))
        #    MP2_check[gpss[m]] = mean(MP2_dict[gpss[m]])
        #MP2_check["TOTAL"] = mean(MP2_values)
        ION_check = {}
        ION_check["TYPE"] = "ION result"
        ION_values1 = []
        ION_values2 = []
        for m in range(len(gpss)):
            amount = len(ION_dict[gpss[m]])
            passed = []
            for n in range(len(ION_dict[gpss[m]])):
                if -10< ION_dict[gpss[m]][n] < 10:
                    passed.append(ION_dict[gpss[m]][n])
                    ION_values1.append(ION_dict[gpss[m]][n])
                ION_values2.append(ION_dict[gpss[m]][n])
            ION_check[gpss[m]] =len(passed)/amount
        ION_check["TOTAL"] = len(ION_values1)/len(ION_values2)
        #ION_values = []
#
        #for m in range(len(gpss)):
        #    for n in range(len(ION_dict[gpss[m]])):
        #        ION_dict[gpss[m]][n] = abs(ION_dict[gpss[m]][n])
        #        ION_values.append(abs(ION_dict[gpss[m]][n]))
        #    ION_check[gpss[m]] = mean(ION_dict[gpss[m]])
        #ION_check["TOTAL"] = mean(ION_values)
        IOD_check = {}
        IOD_check["TYPE"] = "IOD result"
        IOD_values1 = []
        IOD_values2 = []
        for m in range(len(gpss)):
            amount = len(IOD_dict[gpss[m]])
            passed = []
            for n in range(len(IOD_dict[gpss[m]])):
                if -0.3< IOD_dict[gpss[m]][n] < 0.3:
                    passed.append(IOD_dict[gpss[m]][n])
                    IOD_values1.append(IOD_dict[gpss[m]][n])
                IOD_values2.append(IOD_dict[gpss[m]][n])
            IOD_check[gpss[m]] = len(passed) / amount
        IOD_check["TOTAL"] = len(IOD_values1)/len(IOD_values2)
        #IOD_values = []
#
        #for m in range(len(gpss)):
        #    for n in range(len(IOD_dict[gpss[m]])):
        #        IOD_dict[gpss[m]][n] = abs(IOD_dict[gpss[m]][n])
        #        IOD_values.append(abs(IOD_dict[gpss[m]][n]))
        #    IOD_check[gpss[m]] = mean(IOD_dict[gpss[m]])
        #IOD_check["TOTAL"] = mean(IOD_values)
        cyc_check = {}
        cyc_check["TYPE"] = "cyc result"
        cyc_values1 = []
        cyc_values2 = []
        for m in range(len(gpss)):
            amount = len(cyc_dict[gpss[m]])
            passed = []
            for n in range(len(cyc_dict[gpss[m]])):
                if abs(cyc_dict[gpss[m]][n]) < 2:
                    passed.append(cyc_dict[gpss[m]][n])
                    cyc_values1.append(cyc_dict[gpss[m]][n])
                cyc_values2.append(abs(cyc_dict[gpss[m]][n]))
            cyc_check[gpss[m]] = len(passed)/amount
        cyc_check["TOTAL"] = len(cyc_values1)/len(cyc_values2)
        print(ele_check)
        print(MP1_check)
        print(MP2_check)
        print(ION_check)
        print(IOD_check)
        print(cyc_check)
        fieldnames = gpss
        fieldnames.insert(0, "TYPE")
        fieldnames.append("TOTAL")
        row1 = [ele_check, MP1_check, MP2_check, ION_check, IOD_check, cyc_check]

        with open(path[:-4] + "summary.csv", "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(row1)
def ReceiverEditor(path, newpath, receiverlibrary):
    rinex_lines = rt(path)
    rec_type = "".join(rinex_lines[RO(path).rec()][20: 40])
    for k, v in receiverlibrary.items():
        if "{0:<20}".format(k) == rec_type:
            rinex_lines[RO(path).rec()] = rinex_lines[RO(path).rec()]\
                .replace(rec_type, "{0:<20}".format(v))
            try:
                with open(newpath, "w") as fb:
                    for m in range(len(rinex_lines)):
                            fb.write(rinex_lines[m])
            except UnicodeEncodeError:
                with open(newpath, "w", encoding=dc(path)) as fb:
                    for m in range(len(rinex_lines)):
                            fb.write(rinex_lines[m])
def AntennaEditor(path, newpath):#just for antenna lines' rad not for all
    rinex_lines = rt(path)
    ant_rad = "".join(rinex_lines[RO(path).ant()][36: 40])
    if ant_rad == "    ":
        rinex_lines[RO(path).ant()] = rinex_lines[RO(path).ant()]\
                .replace(ant_rad, "NONE")
        try:
            with open(newpath, "w") as fb:
                for m in range(len(rinex_lines)):
                    fb.write(rinex_lines[m])
        except UnicodeEncodeError:
            with open(newpath, "w", encoding=dc(path)) as fb:
                for m in range(len(rinex_lines)):
                    fb.write(rinex_lines[m])
def DataFinding(rootpath, keywords1="", RINEXextension=""):
    all_files = []
    allfiles = sf(rootpath, all_files)
    RINEX_FILES = []
    for n in range(len(allfiles)):
        if keywords1 in allfiles[n] and RINEXextension.upper() in allfiles[n] or \
                keywords1 in allfiles[n] and RINEXextension.lower() in allfiles[n]:
            RINEX_FILES.append((allfiles[n]))
    return RINEX_FILES
def date2doy(year, month, day):
    month_leapyear = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    month_notleap = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    doy = 0
    if month == 1:
        pass
    elif year % 4 == 0 and (year % 100 != 0 or year % 400 == 0):
        for i in range(month - 1):
            doy += month_leapyear[i]
    else:
        for i in range(month - 1):
            doy += month_notleap[i]
    doy += day
    return doy
def DataCleaning(RINEX_FILES,ReceiverLibraryPath, AntennaLibraryPath, newfolder_root, reportpath):
    fieldnames = ["origin path", "new path","version", "marker", "longitude", "latitude", "rec type", "ant type"]
    row = []
    ReceiverLibrary = {}
    reclines = rt(ReceiverLibraryPath)
    for n in range(len(reclines)):
        pause_index = reclines[n].index("!")
        ReceiverLibrary[reclines[n][0:pause_index].replace(',', '')] = reclines[n][pause_index+1: -1].replace(',', '')
    AntennaLibrary = {}
    antlines = rt(AntennaLibraryPath)
    for n in range(len(antlines)):
        pause_index = antlines[n].index("!")
        AntennaLibrary[antlines[n][0:pause_index].replace(',', '')] = antlines[n][pause_index + 1: -1].replace(',', '')
    o_path = []
    for n in range(len(RINEX_FILES)):
        if RINEX_FILES[n][-1] == "o" or RINEX_FILES[n][-1] == "O":
            o_path.append(RINEX_FILES[n])

    for n in range(len(o_path)):
        csv_contents = {}
        path = o_path[n]
        print(path)
        name = (o_path[n][(len(list(o_path[n])) - list(o_path[n])[::-1].index("\\")):])
        try:
            lines = rt(path)
            try:
                #n_lines = rt(n_path[n])
                firsttime = RO(path).starttime()
                try:
                    doy = date2doy(int(firsttime.split()[0]), int(firsttime.split()[1]), int(firsttime.split()[2]))
                    marker_index = RO(path).marker_name()
                    new_marker = lines[marker_index][0:4]
                    if " " in (new_marker):
                        new_marker = new_marker.replace(" ", "a")

                    csv_contents['version'] = RO(path).version()
                    csv_contents["marker"] = "*" + str(new_marker)
                    lon_lat = xyz2BLH(RO(path).xyz()[0], RO(path).xyz()[1], RO(path).xyz()[2], rad=False)
                    csv_contents["longitude"] = lon_lat[0]
                    csv_contents["latitude"] = lon_lat[1]
                    csv_contents["origin path"] = path
                    #for m in range(RO(path).end_of_header()):
                    for m in range(20):
                        neirong = list(lines[m])
                        result = [i for i in neirong if not re.findall("[^\u0000-\u05C0\u2100-\u214F]+", i)]
                        if len(result) != len(neirong):
                            for k in range(len(neirong)):
                                if neirong[k] not in result:
                                    neirong[k] = "-"
                                    print(o_path[n])
                            lines[m] = lines[m].replace(lines[m], "".join(neirong))
                    if len(new_marker) != 4:
                        new_marker += "a"
                    lines[marker_index] = "{0:<60}".format(new_marker) + "MARKER NAME" + "\n"
                    rec_type = "".join(lines[RO(path).rec()][20: 40])
                    for k, v in ReceiverLibrary.items():
                        if "{0:<20}".format(k) == rec_type:
                            lines[RO(path).rec()] = lines[RO(path).rec()].replace(rec_type, "{0:<20}".format(v))
                    csv_contents["rec type"] = "".join(lines[RO(path).rec()][20: 40])
                    ant_type = "".join(lines[RO(path).ant()][20: 40])
                    for k, v in AntennaLibrary.items():
                        if "{0:<20}".format(k) == ant_type:
                            lines[RO(path).rec()] = lines[RO(path).ant()].replace(ant_type, "{0:<20}".format(v))
                    csv_contents["rec type"] = "".join(lines[RO(path).rec()][20: 40])
                    csv_contents["ant type"] = "".join(lines[RO(path).ant()][20: 40])
                    ofile_extension = path[-3:]
                    print(ofile_extension)
                    if ofile_extension[-1] == "o":
                        nfile_extension = ofile_extension.replace("o","n")
                    else:
                        nfile_extension = ofile_extension.replace("O", "N")
                    npath = path.replace(path[-3:], nfile_extension)
                    n_lines = rt(npath)
                    newfolder = newfolder_root + "\\" + ofile_extension[:2] +str(doy)
                    os.makedirs(newfolder, exist_ok=True)
                    o_new_path = newfolder_root + "\\"+ ofile_extension[:2]+ str(doy) + "\\" + "".join([new_marker, str(doy), "0.", ofile_extension])
                    folderfile = fnej(newfolder, "."+ofile_extension)
                    csv_contents["new path"] = o_new_path
                    if o_new_path in folderfile:
                        i = 1
                        pluso_new_path = newfolder_root + "\\"+ ofile_extension[:2] + str(doy) + "\\" + "".join([new_marker, str(doy), "0(" + str(i) + ")." + ofile_extension])
                        csv_contents["new path"] = pluso_new_path
                        while pluso_new_path in folderfile:
                            i += 1
                            pluso_new_path = newfolder_root + "\\"+ ofile_extension[:2] + str(doy) + "\\" + "".join(
                                [new_marker, str(doy), "0(" + str(i) + ")." + ofile_extension])
                        pluso_new_path = newfolder_root + "\\"+ ofile_extension[:2] + str(doy) + "\\" + "".join([new_marker, str(doy), "0(" + str(i) + ")." + ofile_extension])
                        with open(pluso_new_path, "w", encoding="utf-8") as fb:
                            for n in range(len(lines)):
                                fb.write(lines[n])
                        n_new_path = pluso_new_path.replace(ofile_extension, nfile_extension)
                        with open(n_new_path, "w", encoding="utf-8") as fb:
                            for n in range(len(n_lines)):
                                fb.write(n_lines[n])
                    else:
                        with open(o_new_path, "w", encoding="utf-8") as fb:
                            for n in range(len(lines)):
                                fb.write(lines[n])
                        n_new_path = o_new_path.replace(ofile_extension, nfile_extension)
                        with open(n_new_path, "w", encoding="utf-8") as fb:
                            for n in range(len(n_lines)):
                                fb.write(n_lines[n])
                except AttributeError:
                    csv_contents["origin path"] = path
                    csv_contents["new path"] = "AttributeError"
                    #row.append(csv_contents)
            except FileNotFoundError:
                csv_contents["origin path"] = path
                csv_contents["new path"] = "FileNotFoundError"
            #row.append(csv_contents)
        except UnicodeDecodeError:
            csv_contents["origin path"] = path
            csv_contents["new path"] = "UnicodeDecodeError"

        row.append(csv_contents)
    with open(reportpath, "w", encoding="utf-8", newline = "") as f:
        writer = csv.DictWriter(f, fieldnames = fieldnames)
        writer.writeheader()
        writer.writerows(row)
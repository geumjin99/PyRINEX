import os
import re
from PyRINEX.bookworm import detectCode as dc
from PyRINEX.bookworm import readtxt as rt
from PyRINEX.bookworm import file_name_extesion_judgement as fnej
from PyRINEX.bookworm import file_name_extesion_judgement_file as fnejf
from PyRINEX.bookworm import show_files as sf
from PyRINEX.bookworm import readtxt_in_encoding as rte
import itertools
import math
from decimal import Decimal

class RinexO():
    def __init__(self, path):
        #self.name = name
        self.path = path
    def version(self):
        lines = rt(self.path)
        return lines[0][5]
    def interval(self):
        lines = rt(self.path)
        for n in range(len(lines)):
            if "INTERVAL" in lines[n]:
                return float(lines[n].split()[0])
    def rec(self):
        lines = rt(self.path)
        for n in range(len(lines)):
            if "REC # / TYPE / VERS" in lines[n]:
                return n
    def ant(self):
        lines = rt(self.path)
        for n in range(len(lines)):
            if "ANT # / TYPE" in lines[n]:
                return n
    def xyz(self):
        coordinate = []
        lines = rt(self.path)
        for n in range(len(lines)):
            if "APPROX POSITION XYZ" in lines[n]:
                X = float(lines[n].split()[0])
                Y = float(lines[n].split()[1])
                Z = float(lines[n].split()[2])
                coordinate = [X, Y, Z]
                return coordinate
    def marker_name(self):
        lines = rt(self.path)
        for n in range(len(lines)):
            if "MARKER NAME" in lines[n]:
                marker_name = lines[n].split()[0]
                return n
    def marker_number(self):
        lines = rt(self.path)
        for n in range(len(lines)):
            if "MARKER NUMBER" in lines[n]:
                marker_number = lines[n].split()[0]
                return marker_number
    def obstime(self):
        lines = rt(self.path)
        for n in range(len(lines)):
            if "TIME OF FIRST OBS" in lines[n] and "TIME OF LAST OBS" in lines[n+1]:
                observ_time = [lines[n][:43], "-", lines[n+1][:43]]
                return observ_time
    def starttime(self):
        lines = rt(self.path)
        for n in range(len(lines)):
            if "TIME OF FIRST OBS" in lines[n]:
                return lines[n]
    def endtime(self):
        lines = rt(self.path)
        for n in range(len(lines)):
            if "TIME OF LAST OBS" in lines[n]:
                return lines[n]
    def end_of_header(self):
        lines = rt(self.path)
        for n in range(len(lines)):
            if "END OF HEADER" in lines[n]:
                return n
    def obstype(self):
        lines = rt(self.path)
        observation_type = []
        for n in range(len(lines)):
            if lines[0][5] == "2":
                if "# / TYPES OF OBSERV" in lines[n]:
                    for m in range(len("".join(lines[n][6:60]).split())):
                        observation_type.append("".join(lines[n][6:60]).split()[m])
            if lines[0][5] == "3":
                if "SYS / # / OBS TYPES" in lines[n]:
                    sys_obstypes = [lines[n].split()[0]]
                    for m in range(len("".join(lines[n][6:60]).split())):
                        sys_obstypes.append("".join(lines[n][6:60]).split()[m])
                    observation_type.append(sys_obstypes)
        return(observation_type)
    def gpsprn(self):
        lines = rt(self.path)
        gpsprn = []
        for n in range(len(lines)):
            if "PRN / # OF OBS" in lines[n]:
                if "GPS" in lines[0]:
                    gpsprn.append(lines[n].split()[0])
                else:
                    if "G" in lines[n].split()[0]:
                        gpsprn.append(lines[n].split()[0])
        for m in range(len(gpsprn)):
            if "G" not in gpsprn[m]:
                if len(gpsprn[m]) == 1:
                    gpsprn[m] = "G0"+gpsprn[m]
                else:
                    gpsprn[m] = "G"+gpsprn[m]
        return gpsprn
class rinex2_obs():
    def __init__(self, path):
        #self.name = name
        self.path = path
    def epoch_index(self):
        name = os.path.basename(self.path)
        yr = "".join(name[-3:-1])
        lines = rt(self.path)
        e_index = []
        for n in range(RinexO(self.path).end_of_header(), len(lines)):
            if yr == lines[n][1:3]:
                epoch_info = {}
                sat_num = int("".join(lines[n][29:32]))
                epoch_info["time"] = " ".join(lines[n].split()[0:6])
                epoch_info["sat_num"] = sat_num
                #e_index.append(epoch_info)
                epoch_line_amount = math.ceil(sat_num/12)
                sat = ""
                for m in range(n, n+epoch_line_amount):
                    sat += lines[m][32:68].strip("\n")
                sat = re.findall(r'.{3}', sat)
                observation_type = RinexO(self.path).obstype()
                sat_lines_occupied = math.ceil(len(observation_type)/5)
                for m in range(len(sat)):
                    if " " in sat[m]:
                        if sat[m][1] == " ":
                            sat[m] = "G0"+sat[m][-1]
                        else:
                            sat[m] = "G"+sat[m][1]+sat[m][2]
                    observ_per_sat = ""
                    for o in range(1, sat_lines_occupied+1):
                        line_contests = lines[n+epoch_line_amount-1+m*sat_lines_occupied+o][:80].strip("\n")
                        if len(line_contests) != 80:
                            for k in range(80-len(line_contests)):
                                line_contests += " "
                        observ_per_sat += line_contests
                    sat_observ = re.findall(r'.{16}', observ_per_sat)
                    for q in range(len(sat_observ)):
                        if sat_observ[q] == "                ":
                            sat_observ[q] = sat_observ[q].replace("                ", "         0.000  ")
                    observ_dictonary = {}
                    for p in range(len(observation_type)):
                        observ_dictonary[observation_type[p]] = sat_observ[p]
                    epoch_info[sat[m]] = observ_dictonary
                e_index.append(epoch_info)
        return e_index
    def gpslist(self):
        gpslist1 = []

        e_index = rinex2_obs(self.path).epoch_index()
        for n in range(len(e_index)):
            for k, v in e_index[n].items():
                if "G" in k and k not in gpslist1:
                    gpslist1.append(k)
        gpslist1.sort()
        return gpslist1
class rinex3_obs():
    def __init__(self, path):
        #self.name = name
        self.path = path
    def epoch_index(self):
        lines = rt(self.path)
        e_index = []
        observation_type = RinexO(self.path).obstype()
        for n in range(RinexO(self.path).end_of_header(), len(lines)):
            if ">" == lines[n][0]:
                epoch_info = {}
                sat_num = int("".join(lines[n][32:35]))
                epoch_info["time"] = " ".join(lines[n].split()[1:7])
                epoch_info["sat_num"] = sat_num
                observation_type = RinexO(self.path).obstype()
                for m in range(n+1, n+int(sat_num)+1):
                    observ_dictonary = {}
                    observ_per_sat = ""
                    for k in range(len(observation_type)):
                        if lines[m][0] == observation_type[k][0]:

                            obstype_amount = len(observation_type[k])-1

                    j = 1
                    l = 3
                    while j <= obstype_amount:
                        observ_per_sat += "".join(lines[m][l:l+16].strip("\n"))
                        l = l+16
                        j = j+1
                    sat_observ = re.findall(r'.{16}', observ_per_sat)
                    if len(sat_observ) != obstype_amount:
                        for n in range(obstype_amount-len(sat_observ)):
                            sat_observ.append("                ")
                    observ_dictonary = {}
                    for q in range(len(sat_observ)):
                        if sat_observ[q] == "                ":
                            sat_observ[q] = sat_observ[q].replace("                ", "         0.000  ")
                    for p in range(len(observation_type)):
                        if lines[m][0] == observation_type[p][0]:
                            for q in range(len(observation_type[p])-1):
                                observ_dictonary[observation_type[p][q+1]] = sat_observ[q]
                        #可能有必要对单一卫星系统的卫星prn记录中的空格的判断写一个判断式
                        epoch_info[lines[m][0:3]] = observ_dictonary
                e_index.append(epoch_info)
        return e_index
    def gpslist(self):
        gpslist1 = []


        e_index = rinex3_obs(self.path).epoch_index()
        for n in range(len(e_index)):
            for k, v in e_index[n].items():
                if "G" in k and k not in gpslist1:
                    gpslist1.append(k)
        gpslist1.sort()
        return gpslist1
class rinex2_nav():
    def __init__(self, path):
        #self.name = name
        self.path = path
    def IONpar(self):
        lines = rt(self.path)
        IONab = []
        IONdata = ""
        for n in range(len(lines)):
            if "ION ALPHA" in lines[n]:
                for m in range(n, n+2):
                    IONdata += lines[m].strip("\n")[2:50]
                ION = re.findall(r'.{12}', IONdata)
                IONab = ION
        if IONab == []:
            IONab.append("NO ION data.")
        return IONab
    def sat_nav(self):
        Nlist = []
        lines = rt(self.path)
        for n in range(len(lines)):
            if "END OF HEADER" in lines[n]:
                header_index = n
        for m in range(header_index+1, len(lines)):


            if lines[m][1] != " ":
                if lines[m][0] == " ":
                    PRN = "G0" + lines[m][1]
                else:
                    PRN = "G" + lines[m][0:2]
                data = ""
                dic = {}
                sondic = {}
                keys = ["EPOCH", "SV clock bias", "SV clock drift", "SV clock drift rate",
                        "IODE Issue of Data, Ephemeris", "Crs", "Delta n", "M0",
                        "Cuc", "e Eccentricity", "Cus", "sqrt(A)",
                        "Toe Time of Ephemeris", "Cic", "OMEGA0", "CIS",
                        "i0", "Crc", "omega", "OMEGA DOT",
                        "IDOT", "Codes on L2 channel", "GPS Week #", "L2 P data flag",
                        "SV accuracy", "SV health", "TGD", "IODC Issue of Data, Clock",
                        "Transmission time of message", "Fit interval"]
                for o in range(m, m+8):
                    data += lines[o][3:].strip("\n")
                Ndata = re.findall(r'.{19}', data)
                #PRN = lines[m][0:3]
                l = min(len(keys), len(Ndata))
                for k in range(l):
                    sondic[keys[k]] = Ndata[k]
                dic[PRN] = sondic
                Nlist.append(dic)
        return Nlist
class rinex3_nav():
    def __init__(self, path):
        self.path = path
    def IONpar(self):
        lines = rt(self.path)
        IONab = []
        IONdata = ""
        for n in range(len(lines)):
            if "GPSA" in lines[n]:
                for m in range(n, n+2):
                    IONdata += lines[m].strip("\n")[5:53]
                ION = re.findall(r'.{12}', IONdata)
                IONab = ION
        if IONab == []:
            IONab.append("NO ION data.")
        return IONab
    def sat_nav(self):
        Nlist = []
        lines = rt(self.path)
        for n in range(len(lines)):
            if "END OF HEADER" in lines[n]:
                header_index = n
        for m in range(header_index+1, len(lines)):
            if lines[m][0:3] != "   ":
                PRN = lines[m][0:3]

                if PRN[1] == " ":
                    PRN = "".join(["G0", PRN[2]])
                data = ""
                dic = {}
                sondic = {}
                keys = ["EPOCH", "SV clock bias", "SV clock drift", "SV clock drift rate",
                        "IODE Issue of Data, Ephemeris", "Crs", "Delta n", "M0",
                        "Cuc", "e Eccentricity", "Cus", "sqrt(A)",
                        "Toe Time of Ephemeris", "Cic", "OMEGA0", "CIS",
                        "i0", "Crc", "omega", "OMEGA DOT",
                        "IDOT", "Codes on L2 channel", "GPS Week #", "L2 P data flag",
                        "SV accuracy", "SV health", "TGD", "IODC Issue of Data, Clock",
                        "Transmission time of message", "Fit interval", "spare1", "spare2"]
                for o in range(m, m+8):
                    data += lines[o][4:].strip("\n")
                Ndata = re.findall(r'.{19}', data)
                #PRN = lines[m][0:3]
                l = min(len(keys), len(Ndata))
                for k in range(l):
                    sondic[keys[k]] = Ndata[k]
                dic[PRN] = sondic
                Nlist.append(dic)
        return Nlist












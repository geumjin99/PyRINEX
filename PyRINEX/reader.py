import re
import codecs
import os
import chardet
import math
import json

def detectCode(path):
	with open(path, "rb") as file:
		data = file.read(200000)
		dicts = chardet.detect(data)
	return dicts["encoding"]
def show_files(path, all_files):
    file_list = os.listdir(path)
    for file in file_list:
        cur_path = os.path.join(path, file)
        if os.path.isdir(cur_path):
            show_files(cur_path, all_files)
        else:
            all_files.append(cur_path)
    return all_files
def readtxt(path):
    try:
        with open(path, "r") as file_object:
            lines = file_object.readlines()
            return lines
    except UnicodeDecodeError:
        with open(path, "r", encoding=detectCode(path)) as file_object:
            lines = file_object.readlines()
            return lines
def oheader(path):
    HeaderInfo = {}
    contents = readtxt(path)
    version = contents[0][5]
    HeaderInfo["version"] = version
    type = contents[0][40]
    HeaderInfo["type"] = type
    obstypes_linenumber = []
    prns = []
    for n in range(len(contents)):
        if "MARKER NAME" in contents[n]:
            HeaderInfo["MARKER_NAME"] = [n, contents[n].split()[0]]
            HeaderInfo["MARKER_NUMBER"] = [n+1, contents[n+1].split()[0]]
        if "REC # / TYPE / VERS" in contents[n]:
            HeaderInfo["receiver_type"] = [n, contents[n][20:40]]
        if "ANT # / TYPE" in contents[n]:
            HeaderInfo["antenna_type"] = [n, contents[n][20:40]]
        if "APPROX POSITION XYZ" in contents[n]:
            X = float(contents[n].split()[0])
            Y = float(contents[n].split()[1])
            Z = float(contents[n].split()[2])
            HeaderInfo["APPROX POSITION XYZ"] = [X, Y, Z]
        if "TIME OF FIRST OBS" in contents[n]:
            HeaderInfo["TIME OF FIRST OBS"] = ["-".join(contents[n].split()[0:6])]
        if "TIME OF LAST OBS" in contents[n]:
            HeaderInfo["TIME OF FIRST OBS"] = ["-".join(contents[n].split()[0:6])]
        if "# / TYPES OF OBSERV" in contents[n]:
            obstypes_linenumber.append(n)
        if "SYS / # / OBS TYPES" in contents[n]:
            obstypes_linenumber.append(n)
        if "PRN / # OF OBS" in contents[n]:
            prn = contents[n][3:6]
            if prn != "   ":
                prns.append(prn)
        if "END OF HEADER" in contents[n]:
            HeaderInfo["END OF HEADER"] = n
    if prns == []:
        EndOfHeader = HeaderInfo["END OF HEADER"]
        yr = "".join(path[-3:-1])
        if HeaderInfo["version"] == "2":
            for n in range(EndOfHeader, len(contents)):
                if yr == contents[n][1:3] and contents[n][3] == " ":
                    sat_num = int("".join(contents[n][29:32]))
                    epoch_line_amount = math.ceil(sat_num / 12)
                    sat = ""
                    for m in range(n, n + epoch_line_amount):
                        sat += contents[m][32:68].strip("\n")
                    sat = re.findall(r'.{3}', sat)
                    for m in range(len(sat)):
                        if " " in sat[m]:
                            if sat[m][1] == " ":
                                sat[m] = "G0" + sat[m][-1]
                            else:
                                sat[m] = "G" + sat[m][1] + sat[m][2]
                    for m in range(len(sat)):
                        if sat[m] not in prns:
                            prns.append(sat[m])
            prns.sort()
        if HeaderInfo["version"] == "3":
            for n in range(EndOfHeader, len(contents)):
                if contents[n][0] == ">":
                    sat = []
                    i = 1
                    while contents[n+i][0] != ">":
                        sat.append(contents[n+i][:3])
                        #i = i+1
                        if n+i+1 < len(contents):
                            i = i+1
                    for m in range(len(sat)):
                        if sat[m] not in prns:
                            prns.append(sat[m])
            prns.sort()
    for n in range(len(prns)):
        if prns[n][1] == " ":
            prns[n] = "G0" + prns[n][-1]
        if prns[n][0] == " ":
            prns[n] = "G" + prns[n][1:3]
    HeaderInfo["PRNS"] = prns
    for n in range(len(obstypes_linenumber)):
        if "# / TYPES OF OBSERV" in contents[obstypes_linenumber[0]]:
            obslist = []
            for m in range(len(obstypes_linenumber)):
                for o in range(len(contents[obstypes_linenumber[m]][6:60].split())):
                    obslist.append(contents[obstypes_linenumber[m]][6:60].split()[o])
            HeaderInfo["ObsTypes"] = obslist
        if "SYS / # / OBS TYPES" in contents[obstypes_linenumber[0]]:
            obstypes = {}
            for m in range(len(obstypes_linenumber)):
                obslist = []
                for o in range(len(contents[obstypes_linenumber[m]][6:60].split())):
                    obslist.append(contents[obstypes_linenumber[m]][6:60].split()[o])

                obstypes[contents[obstypes_linenumber[m]][0]] = obslist
            HeaderInfo["ObsTypes"] = obstypes
    return json.dumps(HeaderInfo)
def observations(path):
    HeaderInfo = json.loads(oheader(path))
    EndOfHeader = HeaderInfo["END OF HEADER"]
    yr = "".join(path[-3:-1])
    lines = readtxt(path)
    e_index = {}
    if HeaderInfo["version"] == "2":
        observation_type = HeaderInfo["ObsTypes"]
        sat_lines_occupied = math.ceil(len(observation_type) / 5)
        for n in range(EndOfHeader, len(lines)):
            if yr == lines[n][1:3] and lines[n][3] == " ":
                epoch_info = {}
                sat_num = int("".join(lines[n][29:32]))
                epoch_info["sat_num"] = sat_num
                epoch_line_amount = math.ceil(sat_num / 12)
                sat = ""
                for m in range(n, n + epoch_line_amount):
                    sat += lines[m][32:68].strip("\n")

                sat = re.findall(r'.{3}', sat)
                sat = [item for item in sat if item.strip() != ""]

                for m in range(len(sat)):
                    if " " in sat[m]:
                        if sat[m][1] == " ":
                            sat[m] = "G0" + sat[m][-1]
                        else:
                            sat[m] = "G" + sat[m][1] + sat[m][2]
                    observ_per_sat = ""

                    for o in range(1, sat_lines_occupied + 1):
                        line_contests = lines[n + epoch_line_amount - 1 + m * sat_lines_occupied + o][:80].strip("\n")
                        if len(line_contests) != 80:
                            for k in range(80 - len(line_contests)):
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
                e_index[" ".join(lines[n].split()[0:6])]= epoch_info
    if HeaderInfo["version"] == "3":
        observation_type = HeaderInfo["ObsTypes"]
        for n in range(EndOfHeader, len(lines)):
            if lines[n][0] == ">":
                epoch_info = {}
                sat_num = int("".join(lines[n][32:35]))
                epoch_info["sat_num"] = sat_num
                for m in range(n+1, n+1+int(sat_num)):
                    observ_dictonary = {}
                    observ_per_sat = ""
                    PRNcode = lines[m][0]
                    obstypes = observation_type[PRNcode]
                    obstypes_amount = len(obstypes)
                    j = 1
                    l = 3
                    while j <= obstypes_amount:
                        observ_per_sat += "".join(lines[m][l:l + 16].strip("\n"))
                        l = l + 16
                        j = j + 1
                    sat_observ = re.findall(r'.{16}', observ_per_sat)
                    if len(sat_observ) != obstypes_amount:
                        for k in range(obstypes_amount-len(sat_observ)):
                            sat_observ.append("                ")
                    for q in range(len(sat_observ)):
                        if sat_observ[q] == "                ":
                            sat_observ[q] = sat_observ[q].replace("                ", "         0.000  ")
                    for p in range(len(obstypes)):
                        observ_dictonary[obstypes[p]] = sat_observ[p]
                    epoch_info[lines[m][0:3]] = observ_dictonary
                e_index[" ".join(lines[n].split()[1:7])]= epoch_info
    return json.dumps(e_index)
def nheader(path):
    contents = readtxt(path)
    if contents[0][5] == "2":
        IONab = []
        IONdata = ""
        for n in range(len(contents)):
            if "ION ALPHA" in contents[n]:
                for m in range(n, n + 2):
                    IONdata += contents[m].strip("\n")[2:50]
                ION = re.findall(r'.{12}', IONdata)
                IONab = ION
        if IONab == []:
            IONab.append("NO ION data.")
    if contents[0][5] == "3":
        contents = readtxt(path)
        IONab = []
        IONdata = ""
        for n in range(len(contents)):
            if "GPSA" in contents[n]:
                for m in range(n, n + 2):
                    IONdata += contents[m].strip("\n")[5:53]
                ION = re.findall(r'.{12}', IONdata)
                IONab = ION
        if IONab == []:
            IONab.append("NO ION data.")
    return IONab
def navigations(path):
    Nlist = []
    contents = readtxt(path)
    if contents[0][5] == "2":
        for n in range(len(contents)):
            if "END OF HEADER" in contents[n]:
                header_index = n
        for m in range(header_index + 1, len(contents)):
            if contents[m][1] != " ":
                if contents[m][0] == " ":
                    PRN = "G0" + contents[m][1]
                else:
                    PRN = "G" + contents[m][0:2]
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
                for o in range(m, m + 8):
                    data += contents[o][3:].strip("\n")
                Ndata = re.findall(r'.{19}', data)
                #PRN = lines[m][0:3]
                l = min(len(keys), len(Ndata))
                for k in range(l):
                    sondic[keys[k]] = Ndata[k]
                dic[PRN] = sondic
                Nlist.append(dic)
    if contents[0][5] == "3":
        for n in range(len(contents)):
            if "END OF HEADER" in contents[n]:
                header_index = n
        for m in range(header_index+1, len(contents)):
            if contents[m][0:3] != "   ":
                PRN = contents[m][0:3]
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
                    data += contents[o][4:].strip("\n")
                Ndata = re.findall(r'.{19}', data)
                #PRN = contents[m][0:3]
                l = min(len(keys), len(Ndata))
                for k in range(l):
                    sondic[keys[k]] = Ndata[k]
                dic[PRN] = sondic
                Nlist.append(dic)
    dict = {}
    prns = []
    for n in range(len(Nlist)):
        for k,v in Nlist[n].items():
            if k not in prns:
                prns.append(k)
    for n in range(len(prns)):
        dict[prns[n]] = []
        for m in range(len(Nlist)):
            for k,v in Nlist[m].items():
                if k == prns[n]:
                    dict[prns[n]].append(Nlist[m][k])
    return json.dumps(dict)
















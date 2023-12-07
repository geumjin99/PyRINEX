import re
import codecs
import os
import chardet
import csv
from reader import *
import numpy as np


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
def readtxt_in_encoding(path, ec):
    with open(path, "r", encoding=ec) as file_object:
        lines = file_object.readlines()
        return lines
def file_name_extesion_judgement(file_dir, extension_name):
    file_list_withdir = []
    for root, dirs, files in os.walk(file_dir):
        for file in files:
            if os.path.splitext(file)[1] == extension_name.upper() or \
                    os.path.splitext(file)[1] == extension_name.lower():
                file_list_withdir.append(os.path.join(root, file))
        return file_list_withdir
def file_name_extesion_judgement_file(file_dir, extension_name):
    file_list = []
    for root, dirs, files in os.walk(file_dir):
        for file in files:
            if os.path.splitext(file)[1] == extension_name.upper() or \
                    os.path.splitext(file)[1] == extension_name.lower():
                file_list.append(file)
        return file_list
def DataFinding(root_path, keywordslist, RINEXextension):
    RINEXfiles = []
    for foldername, subfolders, filenames in os.walk(root_path):
        for filename in filenames:
            if all(keyword in filename for keyword in keywordslist) and filename.endswith(RINEXextension):
                file_path = os.path.join(foldername, filename)
                RINEXfiles.append(file_path)
    return RINEXfiles
def generate_unique_filename(path):
    if not os.path.exists(path):
        return path
    filename, extension = os.path.splitext(path)
    counter = 1
    while os.path.exists(path):
        path = f"{filename}_{counter}{extension}"
        counter += 1
    return path
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
def DataCleaning(RINEX_FILES,ReceiverLibraryPath, AntennaLibraryPath, newfolder_root):
    fieldnames = ["origin path", "version", "non English", "origin marker", "origin rec", "origin ant",
                  "new path", "marker", "longitude", "latitude", "rec type", "ant type"]
    row = []
    ReceiverLibrary = {}
    AntennaLibrary = {}
    reclines = readtxt(ReceiverLibraryPath)
    for n in range(len(reclines)):
        pause_index = reclines[n].index("?")
        ReceiverLibrary[reclines[n][0:pause_index].replace(',', '')] = reclines[n][pause_index+1: -1].replace(',', '')
    antlines = readtxt(AntennaLibraryPath)
    for n in range(len(antlines)):
        pause_index = antlines[n].index("?")
        AntennaLibrary[antlines[n][0:pause_index].replace(',', '')] = antlines[n][pause_index+1: -1].replace(',', '')
    o_path = RINEX_FILES
    n_path = []
    for n in range(len(o_path)):
        if o_path[n][-1] == "o":
            n_path.append(o_path[n][:-1]+"n")
        else:
            n_path.append(o_path[n][:-1]+"N")

    for n in range(len(o_path)):
        try:
            csv_contents = {}
            path = o_path[n]
            lines = readtxt(path)
            header = json.loads(oheader(path))
            name = os.path.basename(path)
            csv_contents["version"] = header["version"]
            firsttime = header["TIME OF FIRST OBS"][0].replace("-", " ")
            try:
                doy = date2doy(int(firsttime.split()[0]), int(firsttime.split()[1]), int(firsttime.split()[2]))
            except AttributeError:
                doy = "xxx"
            marker_index = header["MARKER_NAME"][0]
            try:
                csv_contents["origin marker"] = header["MARKER_NAME"][1]
            except TypeError:
                csv_contents["origin marker"] = "x"
            new_marker = header["MARKER_NAME"][1][0:4]
            if " " in new_marker:
                new_marker.replace(" ", "a")
            if len(new_marker) != 4:
                while len(new_marker) < 4:
                    new_marker += "a"
            csv_contents["marker"] = str(new_marker)
            lines[marker_index] = "{0:<60}".format(new_marker) + "MARKER NAME" + "\n"
            x = header["APPROX POSITION XYZ"][0]
            y = header["APPROX POSITION XYZ"][1]
            z = header["APPROX POSITION XYZ"][2]
            lon_lat = xyz2BLH(x, y, z, rad=False)
            csv_contents["longitude"] = lon_lat[0]
            csv_contents["latitude"] = lon_lat[1]
            csv_contents["origin path"] = path
            for m in range(20):
                neirong = list(lines[m])
                result = [i for i in neirong if not re.findall("[^\u0000-\u05C0\u2100-\u214F]+", i)]
                if len(result) != len(neirong):
                    csv_contents["non English"] = "Non-English characters present"
                    for k in range(len(neirong)):
                        if neirong[k] not in result:
                            neirong[k] = "-"
                    lines[m] = lines[m].replace(lines[m], "".join(neirong))
            rec_type = header["receiver_type"][1]
            csv_contents["origin rec"] = lines[header["receiver_type"][0]][20:40]
            for k, v in ReceiverLibrary.items():
                if "{0:<20}".format(k) == rec_type:
                    lines[header["receiver_type"][0]] = lines[header["receiver_type"][0]].replace(rec_type, "{0:<20}".format(v))
            csv_contents["rec type"] = lines[header["receiver_type"][0]][20:40]
            ant_type = header["antenna_type"][1]
            for k, v in AntennaLibrary.items():
                if "{0:<20}".format(k) == ant_type:
                    lines[header["antenna_type"][0]] = lines[header["antenna_type"][0]].replace(ant_type, "{0:<20}".format(v))
            csv_contents["ant type"] = lines[header["antenna_type"][0]][20:40]
            ant_rad = header["antenna_type"][1][36:40]
            csv_contents["origin ant"] = lines[header["antenna_type"][0]][20:40]
            if len(ant_rad) == 0:
                lines[header["antenna_type"][0]] = "".join(
                    [lines[header["antenna_type"][0]][:36], "NONE", lines[header["antenna_type"][0]][40:]])
            csv_contents["ant type"] = lines[header["antenna_type"][0]][20:40]
            newfolder = newfolder_root + "\\" + path[-3:-1] + str(doy)
            os.makedirs(newfolder, exist_ok=True)
            o_new_path = newfolder_root + "\\" + path[-3:-1] + str(doy) + "\\" + "".join([new_marker, str(doy), "0.", path[-3:]])
            o_new_path = generate_unique_filename(o_new_path)
            if o_new_path[-1] == "o":
                n_new_path = o_new_path[:-1] + "n"
            else:
                n_new_path = o_new_path[:-1] + "N"
            csv_contents["new path"] = o_new_path
            with open(o_new_path, "w", encoding="utf-8") as fb:
                for n1 in range(len(lines)):
                    fb.write(lines[n1])
            status = os.path.exists(n_path[n])
            if status != False:
                nlines = readtxt(n_path[n])
                print(n_path[n])
                with open(n_new_path, "w", encoding="utf-8") as fb:
                    for n2 in range(len(nlines)):
                        fb.write(nlines[n2])
        except UnicodeDecodeError:
            csv_contents["origin path"] = path
            csv_contents["new path"] = "UnicodeDecodeError"
        except KeyError:
            csv_contents["origin path"] = path
            csv_contents["new path"] = "KeyError"
        row.append(csv_contents)
    with open(newfolder_root + "\\" +"report.csv", "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames = fieldnames)
        writer.writeheader()
        writer.writerows(row)



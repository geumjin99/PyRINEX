import re
import codecs
import os
import chardet

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


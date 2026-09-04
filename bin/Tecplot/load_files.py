#!/usr/bin/env python

import tecplot as tp
from tecplot.exception import *
from tecplot.constant import *
import easygui
import sys
import glob
import os


script_directory = os.path.dirname(os.path.abspath(__file__))
directory = easygui.diropenbox(msg="Locate Data Dir", title=None, default=os.getcwd())
#myglob = easygui.enterbox("Enter glob string",default="data*.dat")
#easygui.msgbox(directory, 'Directory Name')

# Add variable names for Tecplot if they are not already present.
variables_header = 'VARIABLES = "x" "y" "u" "v" "T" "p"'
files = glob.glob(os.path.join(directory, "data*_b0.dat"))
for file_path in files:
    with open(file_path, encoding="utf-8") as data_file:
        contents = data_file.read()
    if not contents.startswith(variables_header):
        contents = variables_header + "\n" + contents
    with open(file_path, "w", encoding="utf-8") as data_file:
        data_file.write(contents)

# Add variable names for Tecplot if they are not already present.
variables_header = 'VARIABLES = "x" "y" "T"'
files = glob.glob(os.path.join(directory, "data*_b1.dat"))
for file_path in files:
    with open(file_path, encoding="utf-8") as data_file:
        contents = data_file.read()
    if not contents.startswith(variables_header):
        contents = variables_header + "\n" + contents
    with open(file_path, "w", encoding="utf-8") as data_file:
        data_file.write(contents)

# Uncomment the following line to connect to a running instance of Tecplot 360:
tp.session.connect()

tp.load_layout(os.path.join(script_directory, "Tecplot", "solid-liquid.lay"))
bnum = 0
ReadDataOption = tp.constant.ReadDataOption.Replace
while True:
    globstring = f"data*_b{bnum}.dat"
    files = glob.glob(directory + "/" + globstring)
    if len(files) == 0:
        break
    tp.data.load_tecplot(files,read_data_option=ReadDataOption)
    ReadDataOption = tp.constant.ReadDataOption.Append
    bnum += 1

frame = tp.active_frame()
frame.load_stylesheet(os.path.join(script_directory, "Tecplot", "solid-liquid.sty"))

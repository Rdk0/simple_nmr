""" https://www.mfitzp.com/article/1d-1h-nmr-data-processing/#loading%20bruker%20fid%20spectra
https://github.com/jjhelmus/nmrglue/tree/master/examples/bruker_data
A basic program plotting processed Bruker data (after ff, phase correction and reference stting) using nmrglue and
matplotlib.  Does not have an integration functionality
"""
import matplotlib.pyplot as plt
import nmrglue as ng
import numpy as np
from tkinter import Tk
from tkinter.filedialog import askdirectory
import pandas as pd
import pandas_read_xml as pdx
import os
from typing import Dict, Tuple

plt.style.use("ggplot")


def read_data(full_path_to_pdata_folder: str) -> Tuple[Dict, np.ndarray]:
    """reading Bruker's processed data"""
    dic, data = ng.bruker.read_pdata(full_path_to_pdata_folder)
    return dic, data


def plotspectra(ppms, data, start=None, stop=None) -> None:
    if start:
        ixs = list(ppms).index(start)
        ppms = ppms[ixs:]
        data = data[ixs:]
    if stop:
        ixs = list(ppms).index(stop)
        ppms = ppms[:ixs]
        data = data[:ixs]
    fig = plt.figure(figsize=(12, 4))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(ppms, data)
    ax.set_xlabel("ppm")
    ax.set_title("1D NMR Spectrum")
    ax.invert_xaxis()


def calculate_ppms(dic: Dict) -> np.ndarray:
    """calculate chemical shifts"""
    zero_fill_size = 65536
    offset = (float(dic["acqus"]["SW"]) / 2) - (
        float(dic["acqus"]["O1"]) / float(dic["acqus"]["BF1"])
    )
    start = float(dic["acqus"]["SW"]) - offset
    end = -offset
    step = float(dic["acqus"]["SW"]) / zero_fill_size
    ppms = np.arange(start, end, -step)[:zero_fill_size]
    return ppms


def get_pdata_folder_name(parent_folder: str) -> str:
    """
    use TKinter to interactively retrieve a folder full path to the processed files: typically ending in  \pdata\1
    """
    title: str = "find the pdata folder"
    Tk().withdraw()
    folder_name: str = askdirectory(initialdir=parent_folder, title=title)
    return folder_name


def read_peak_list(
    peaklist_xml: str, threshold_percentile: float = 0.9
) -> pd.DataFrame:
    """
    returns peaks above certain percentile
    :param peaklist_xml: file name
    :return: data frame
    """
    df: pd.DataFrame = pdx.read_xml(peaklist_xml, ["PeakList", "PeakList1D", "Peak1D"])
    df = df.astype({"@F1": "float", "@intensity": float})
    percentile: float = df["@intensity"].quantile(threshold_percentile)
    df = df.rename(columns={"@F1": "ppm", "@intensity": "intensity"})
    return df[df["intensity"].gt(percentile)]


def main():
    parent_folder: str = r"C:\instruments\data2021\nmr"
    full_path_to_pdata_folder: str = get_pdata_folder_name(parent_folder)
    dic, data = read_data(full_path_to_pdata_folder)
    ppms: np.darray = calculate_ppms(dic)
    peaklist_xml_file: str = os.path.join(full_path_to_pdata_folder, "peaklist.xml")
    df_peaks: pd.DataFrame = read_peak_list(
        peaklist_xml=peaklist_xml_file, threshold_percentile=0.9
    )
    print(df_peaks[["ppm", "intensity"]])
    plotspectra(ppms, data)


if __name__ == "__main__":
    main()

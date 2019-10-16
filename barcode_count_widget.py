#!/usr/bin/env python3
"""
BARET
BARcode Extraction Tool

Author: Jack Feldman
Fall 2019
"""

import subprocess
import os
from tkinter import *
from tkinter import messagebox, scrolledtext
from tkinter.filedialog import askopenfilename, askdirectory
import pandas as pd
import datetime

from modules.extract import *
from modules.metrics import *

# -------------------------------------------------------------------
# BUILD WINDOW
# -------------------------------------------------------------------
window = Tk()
window.title('BARET (BARcode Extraction Tool)')
window.geometry('800x610')
window.configure(background='deepskyblue2')

# FASTQ FILE DIRECTORY PATH
text = Label(window, text='Path to directory containing fastq.gz files:', background='deepskyblue2')
text.grid(row=0, column=0, sticky=W)
fastq = Entry(window, width=30, highlightbackground='deepskyblue2')
fastq.grid(row=0, column=1)

def get_fastq_file():
    filename = askdirectory()
    fastq.delete(0, END)
    fastq.insert(END, filename)

fastq_button = Button(window, text='Select', command=get_fastq_file)
fastq_button.grid(row=0, column=2)


# BARCODE FILE AND SHEET NAME
text = Label(window, text='Path to barcode csv file:', background='deepskyblue2')
text.grid(row=2, column=0, sticky=W)
barcode = Entry(window, width=30, highlightbackground='deepskyblue2')
barcode.insert(END, 'barcodes.csv')
barcode.grid(row=2, column=1)

def get_barcode_file():
    filename = askopenfilename()
    barcode.delete(0, END)
    barcode.insert(END, filename)

barcode_button = Button(window, text='Select', command=get_barcode_file)
barcode_button.grid(row=2, column=2)

# OUTPUT DIRECTORY NAME
text = Label(window, text="Output directory (will be created if doesn't exist):", background='deepskyblue2')
text.grid(row=1, column=0, sticky=W)
output_ = Entry(window, width=30, highlightbackground='deepskyblue2')
output_.grid(row=1, column=1)

def get_output_dir():
    filename = askdirectory()
    output_.delete(0, END   )
    output_.insert(END, filename)

output_button = Button(window, text='Select', command=get_output_dir)
output_button.grid(row=1, column=2)

# OUTPUT TEXT
output_text = scrolledtext.ScrolledText(window,width=50,height=30)
output_text.grid(row=7, column=1, pady=20,)

# -------------------------------------------------------------------
# END BUILD WINDOW
# -------------------------------------------------------------------

def check_inputs(fastq,
                 barcode,
                 output
                 ) -> bool:
    # FASTQ CHECKS
    if fastq == '':
        messagebox.showinfo('Error', 'You must enter a file path.')
    elif not os.path.isdir(fastq):
        messagebox.showinfo('Error', 'FASTQ directory not found.')
    # BARCODE FILE CHECKS
    elif barcode == '':
        messagebox.showinfo('Error', 'You must enter a barcode file.')
    elif not os.path.isfile(barcode):
        messagebox.showinfo('Error', 'Barcode file not found.')
    # OUTPUT DIR CHECKS
    elif output == '':
        messagebox.showinfo('Error', 'You must enter an output directory.')
    else:
        return True


# WHAT HAPPENS WHEN YOU CLICK THE BUTTON
def clicked():
    start = time.time()
    output_text.delete(1.0, 1.0)
    fastq_dir = fastq.get()
    barcode_file = barcode.get()
    output_dir = output_.get()

    inputs_ok = check_inputs(
                    fastq = fastq_dir,
                    barcode = barcode_file,
                    output = output_dir,
    )

    if inputs_ok:

        try:
            reads = get_read_paths(fastq_dir)
            barcodes = get_barcodes(barcode_file)
            data = pd.DataFrame(index = barcodes)
            data.index.name = 'barcodes'
            log_data = pd.DataFrame(index = ["Number of instances with 'N' in barcode region",
                                             "Number of instances where barcode not found in barcode region"])
            log_data.index.name = 'Type of error'
            pcr_df = return_pcr_bias_df()

            file_count = len(reads)
            index = 1

            output_text.insert(1.0, "Starting analysis...\n")
            output_text.update()

            # PROCESS FILES
            for read in reads:
                try:
                    col_name = read.split('/')[-1].split('_')[0]
                    output_text.insert(1.0, f"Processing {col_name} ({index} of {file_count})\n")
                    output_text.update()

                    barcode_counts, pcr_counts, error_log = calc_bc_count(read, barcodes)
                    data[col_name] = barcode_counts.values()
                    log_data[col_name] = error_log.values()
                    for base in pcr_counts:
                        pcr_df[base] += pcr_counts[base]

                    index += 1

                except Exception as e:
                    output_text.insert(1.0, str(e) + '\n')
                    output_text.update()
                    return 1

            output_text.insert(1.0, f"File parsing complete.\n")
            output_text.update()
            data.reset_index(inplace=True)

            # CREATE OUTPUT DIR IF DOESN'T EXIST
            if os.path.isdir(output_dir):
            	pass
            else:
                output_text.insert(1.0, f"Creating directory {output_dir}\n")
                output_text.update()
                os.mkdir(output_dir)

            # WRITE RAW COUNTS AND LOG FILE
            raw_counts_path = os.path.join(output_dir, 'rawcounts.csv')
            norm_counts_path = os.path.join(output_dir, 'normcounts.csv')
            log_name = "logfile_" + datetime.datetime.now().strftime("%Y%m%d") + ".log"
            log_path = os.path.join(output_dir, log_name)

            data.to_csv(raw_counts_path)
            log_data.to_csv(log_path)

            for i in range(11):
                pcr_df.iloc[i] = round(pcr_df.iloc[i] / sum(pcr_df.iloc[i]), 2)
            pcr_path = os.path.join(output_dir, 'pcr_bias.csv')
            pcr_df.to_csv(pcr_path)

            # RUN R NORMALIZATION SCRIPT
            output_text.insert(1.0, f"Normalizing counts...\n")
            output_text.insert(1.0, '\n')
            output_text.update()

            norm_command = ['Rscript', '--vanilla', './modules/normalization.R', raw_counts_path, output_dir]
            process = subprocess.Popen(norm_command, stdout=subprocess.PIPE)
            output, error = process.communicate()

            if process.returncode != 0:
                output_text.insert(1.0, f"*Error in normalization script*\n")
                output_text.insert(1.0, f"Raw counts file written to:\n{raw_counts_path}\n")
                output_text.update()
                return 1

            else:
                output_text.insert(1.0, output)
                output_text.insert(1.0, '\n')
                output_text.update()

            # Calculate avg and std of norm counts for all cell types
            cell_type_variance, normalized_compiled = calc_metrics(data=norm_counts_path)
            variance_path = os.path.join(output_dir, 'cell_type_variance.csv')
            cell_type_variance.to_csv(variance_path)

            normalized_compiled_path = os.path.join(output_dir, 'normalized_compiled.csv')
            normalized_compiled.to_csv(normalized_compiled_path)

            output_text.insert(1.0, f"Files written to:\n{output_dir}\n")
            end = time.time()
            output_text.insert(1.0, "Total Runtime: "+str(datetime.timedelta(seconds=end-start))+"\nAnalysis complete.\n")

        except Exception as e:
            output_text.insert(1.0, "Error: " + str(e) + '\n')
            return 1

        else:
            return 0

start_button = Button(window, text='Start', command=clicked, height = 3, width=20)
start_button.grid(row=7, column=0, sticky=N, pady=20, )

if __name__ == "__main__":
    window.mainloop()

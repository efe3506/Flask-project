"""
This script runs the application using a development server.
It contains the definition of routes and views for the application.
"""

from flask import Flask, render_template, request, redirect, url_for
from flask import send_file
from flask import current_app
from flask import send_from_directory
app = Flask(__name__)



@app.route('/')
def index():
    return render_template('mainPage.html')

def difexp_TF(file_name, out_name):     
    import pandas as pd 
    import math as mt
    folder = pd.read_csv(file_name, sep="\t")

    X = "TF"      # Daha evrensel bir çözüm bulunmalı. Direkt dosyalarda bu sütun X olarak isimlendirilebilir. 

    for i, j in folder.iterrows():
        folder.loc[i, "log(FC)"] = mt.log((j["Tumor"]/j["Normal"]), 2)
    thrs = 1.5
    folder_DE = folder[abs(folder["log(FC)"])> thrs]
    folder_DE_sorted = folder_DE.sort_values("log(FC)", axis=0, ascending= False).loc[:, [X, "log(FC)"]] 
    for k, l in folder_DE_sorted.iterrows():
        if l["log(FC)"]>0:
            folder_DE_sorted.loc[k, "label"] = "+"    
        else: 
            folder_DE_sorted.loc[k, "label"] = "-"

            
    folder_DE_sorted.to_csv(out_name, sep="\t")
    return folder_DE_sorted
    
def difexp_miRNA(file_name, out_name):     # log ve label degerlerinin miRNA icin alinmasi
    import pandas as pd 
    import math as mt
    folder = pd.read_csv(file_name, sep="\t")

    X = "miRNA"      # Daha evrensel bir çözüm bulunmalı. Direkt dosyalarda bu sütun X olarak isimlendirilebilir. 

    for i, j in folder.iterrows():
        folder.loc[i, "log(FC)"] = mt.log((j["Tumor"]/j["Normal"]), 2)
    thrs = 1.5
    folder_DE = folder[abs(folder["log(FC)"])> thrs]
    folder_DE_sorted = folder_DE.sort_values("log(FC)", axis=0, ascending= False).loc[:, [X, "log(FC)"]] 
    for k, l in folder_DE_sorted.iterrows():
        if l["log(FC)"]>0:
            folder_DE_sorted.loc[k, "label"] = "+"    
        else: 
            folder_DE_sorted.loc[k, "label"] = "-"
    folder_DE_sorted.to_csv(out_name, sep="\t")
    return folder_DE_sorted

def mergeall(file1, file2, file3, file4, file5):            #TF-miRNA sonuclarinin merge edilme kismi
    TF_DE_sorted = difexp_TF(file1, file2)
    miRNA_DE_sorted =difexp_miRNA(file3, file4)
    f=open(file5, "w")
    f.write("TF\tup_down\tmiRNA\tup_down\n")

    for i, j in TF_DE_sorted.iterrows():
        for k, l in miRNA_DE_sorted.iterrows():
            if (j["label"] == "+") & (l["label"] == "+") | (j["label"] == "-") & (l["label"] == "-") :
                f.write(j["TF"] + "\t" + j["label"] + "\t" + str(l["miRNA"]) + "\t" + l["label"] + "\n")
    f.close()

def intersect1(input_1, input_2, output):           #ikililer icin intersect
    import pandas as pd 
    sorted_1 = (pd.read_csv(input_1, sep = "\t")).sort_values("TF")
    sorted_2 = (pd.read_csv(input_2, sep = "\t")).loc[:, ["TF", "miRNA"]].sort_values("TF")
    intersect = (pd.merge(sorted_1, sorted_2, how = "inner")).drop_duplicates(keep = "first")
    

    TF_miRNA = pd.read_csv(output, sep = "\t")
    TF_mRNA = pd.read_csv("TF_mRNA_ornek.csv", sep = "\t")
    intersect_2 = pd.merge(TF_miRNA, TF_mRNA, how = "inner")
    intersect_2  

    intersect.to_csv(output, sep = "\t")
    return intersect_2
    
def intersect2(input_1, input_2, output, file6):            #ucluler icin intersect
    import pandas as pd 
    sorted_1 = (pd.read_csv(input_1, sep = "\t")).sort_values("TF")
    sorted_2 = (pd.read_csv(input_2, sep = "\t")).loc[:, ["TF", "miRNA"]].sort_values("TF")
    intersect = (pd.merge(sorted_1, sorted_2, how = "inner")).drop_duplicates(keep = "first")
    

    TF_miRNA = pd.read_csv(output, sep = "\t")
    TF_mRNA = pd.read_csv("TF_mRNA_ornek.csv", sep = "\t")
    intersect_2 = pd.merge(TF_miRNA, TF_mRNA, how = "inner")

    miRNA_mRNA = pd.read_csv("miRNA_mRNA_ornek.csv", sep = "\t")
    intersect_3 = pd.merge(intersect_2, miRNA_mRNA, how = "inner") 
    intersect_3.to_csv("intersect_3.csv", sep = "\t")

    miRNA_TF = pd.read_csv("miRNA_TF_ornek.csv", sep = "\t")
    intersect_4 = pd.merge(intersect_3, miRNA_TF, how = "inner")
    intersect_4.to_csv("intersect_4.csv", sep = "\t")

    f = open(file6, "w")
    f.write("CIRCUITS\n\n\nmRNA\t\t TF\t\tmiRNA")
    for i, j in intersect_4.iterrows():
        f.write("\n\n\n" + "\t\t  |----------\n" + j["mRNA"] + "\t"  + "<---" + "\t" + j["TF"] + "\t" + "--->" + "\t" + j["miRNA"] + "\n"  + "\t" + "|-------------------" + "\n\n\n")
    f.close()

@app.route('/', methods=['POST'])   #file post edildiginde burasi cagirilsin diye
def upload_file():
    TF_expression_file = request.files['file1']
    miRNA_expression_file = request.files['file2']

    TF_expression_file.save(TF_expression_file.filename)
    miRNA_expression_file.save(miRNA_expression_file.filename)

    mergeall(TF_expression_file.filename, "out_1.csv", miRNA_expression_file.filename, "out_2.csv", "merge_all.csv") #TF-miRNA merge etme kismi

    intersect1("TF_miRNA_raw_degistirilmis.csv", "merge_all.csv", "intersect_TF_miRNA.csv") # ikililer icin intersect

    intersect2("TF_miRNA_raw_degistirilmis.csv", "merge_all.csv", "intersect_TF_miRNA.csv", "circuits.txt")   #ucluler icin intersect
                
    #app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0   #to change file content in every file changing       
    return redirect(url_for('index'))  #kendi uzerine donsun diye index yazilir

@app.route('/uploads/<path:filename>', methods=['GET', 'POST']) #to see content of file without downloading
def downloads(filename):
    import csv
    try:
          uploads = os.path.join(current_app.root_path )

          with open('result.txt', "w") as my_output_file:
              with open(filename, "r") as my_input_file:
                  [ my_output_file.write(" ".join(row)+'\n') for row in csv.reader(my_input_file)]
              my_output_file.close()           
          return send_from_directory( directory=uploads, filename = 'result.txt')
    except:
         return redirect(url_for('index')) 

@app.route('/download')
def download_file():        #to download TF file
    try:
        path = "out_1.csv"
        return send_file(path, as_attachment=True)
    except:
        return redirect(url_for('index')) 
@app.route('/download2')
def download_file2():        #to download miRNA file
    try:
        path = "out_2.csv"
        return send_file(path, as_attachment=True)
    except:
        return redirect(url_for('index')) 

if __name__ == '__main__':  # projenin baslatilmasi, portun acilmz
    import os
    HOST = os.environ.get('SERVER_HOST', 'localhost')
    try:
        PORT = int(os.environ.get('SERVER_PORT', '5555'))
    except ValueError:
        PORT = 5555
    app.run(HOST, PORT)


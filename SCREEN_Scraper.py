#Importing libraries
import os
os.chdir("C:/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/")
import requests
import re
import csv
import json
import time
from random import seed
from random import randint
import pandas as pd
df = pd.read_csv('DMR.table.csv')
df
Annotation = df.Annotation
chromosome = df.chr
start = df.start
start[counter]
end = df.end
enh_df = pd.DataFrame([], columns=['chrom', 'start', 'len', 'ctcf_zscore', 'dnase_zscore', 'enhancer_zscore', 'promoter_zscore', 'accession', 'isporximal', 'concordant', 'ctcfmax', 'k4me3max', 'k27acmax', 'vistaids'])

url = 'https://screen-beta-api.wenglab.org/dataws/cre_table'

headers = {
  'authority': 'screen-beta-api.wenglab.org',
  'accept': 'application/json',
  'accept-language': 'en-US,en;q=0.9',
  'content-type': 'application/json',
  'origin': 'https://screen.wenglab.org',
  'referer': 'https://screen.wenglab.org/',
  'sec-ch-ua': '\".Not/A)Brand\";v=\"99\", \"Google Chrome\";v=\"103\", \"Chromium\";v=\"103\"',
  'sec-ch-ua-mobile': '?1',
  'sec-ch-ua-platform': '\"Android\"',
  'sec-fetch-dest': 'empty',
  'sec-fetch-mode': 'cors',
  'sec-fetch-site': 'same-site',
  'user-agent': 'Mozilla/5.0 (Linux; Android 6.0; Nexus 5 Build/MRA58N) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/103.0.0.0 Mobile Safari/537.36'
}
counter = 0
for count in chromosome: 
    data = {
        "uuid": "bd69d494-cc03-45fb-b368-1e7a7c552a28",
        "assembly":"mm10",
        "accessions":[],
        "coord_chrom":str(chromosome[counter]),
        "coord_start":int(start[counter]),
        "coord_end":int(end[counter]),
        "gene_all_start":0,
        "gene_all_end":5000000,
        "gene_pc_start":0,
        "gene_pc_end":5000000,
        "rank_dnase_start":1.64,
        "rank_dnase_end":10,
        "rank_promoter_start":-10,
        "rank_promoter_end":10,
        "rank_enhancer_start":-10,
        "rank_enhancer_end":10,
        "rank_ctcf_start":-10,
        "rank_ctcf_end":10,
        "cellType": None,
        "element_type": None
    }
    url_output = requests.post(url, headers=headers, data=json.dumps(data))
    enh_json = url_output.json()
    enh_list = enh_json.get('cres')
    for idlength, lenth in enumerate(enh_list):
        enh_list_two = enh_list[idlength]['info']
        df_data = [{'chrom':enh_list[idlength]['chrom'], 'start':enh_list[idlength]['start'], 'len':enh_list[idlength]['len'], 'ctcf_zscore':enh_list[idlength]['ctcf_zscore'], 'dnase_zscore':enh_list[idlength]['dnase_zscore'], 'enhancer_zscore':enh_list[idlength]['enhancer_zscore'], 'promoter_zscore':enh_list[idlength]['promoter_zscore'], 'accession':enh_list_two['accession'], 'isproximal':enh_list_two['isproximal'], 'concordant':enh_list_two['concordant'], 'ctcfmax':enh_list_two['ctcfmax'], 'k4me3max':enh_list_two['k4me3max'], 'k27acmax':enh_list_two['k27acmax'], 'vistaids':enh_list[idlength]['vistaids'], 'annotation':Annotation[counter]}]
        enh_df_two = pd.DataFrame(df_data)
        frames = [enh_df, enh_df_two]
        enh_df = pd.concat(frames)
    counter += 1
            
    
enh_df.to_csv('enhancer.csv', index = False)

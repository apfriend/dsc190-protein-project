'''
Written by: Alexander Friend
Last Updated: 2021.05.21

Script to get a .txt file of PDB ids and a 
.json file with the corresponding protein name for each PDB id
'''

import os
import sys
import re
import json
from urllib.request import urlopen, urlretrieve, HTTPError
from bs4 import BeautifulSoup
from tqdm import tqdm
from multiprocessing import Pool

N_PROCESSES=4
PARAMS_FP='../config.json'

def get_html(url):
    '''
    Get html code for webpage from <url>
    ----------
    Parameters
    ----------
        url - str
            url of webpage to get html from
    -------
    Returns
    -------
        html from <url> as a string
    '''
    request=urlopen(url)
    return BeautifulSoup(str(request.read()), features="html.parser")

def find_links(url):
    '''
    Find links to protein entry from <url>
    '''
    #get html of protein page
    html=get_html(url)

    #find all pdb page links
    spans=html.find_all("span",{"class":"rcsb_id_tag"})
    pdb_entries=[]
    for span in spans:
        try:
            href=span.find_all("a")[0]["href"].strip("\\\'")
            pdb_entries.append(href)
        except IndexError as e:
            pass
    #for each pdb page get the pdb file download link
    pdb_links=[]
    for entry in pdb_entries:
        pdb_html=get_html(entry)
        files=pdb_html.find_all("div", id="DownloadFilesButton")[0]
        try:
            download_link=[f["href"] for f in files.find_all("a") if re.match("^.*\.pdb$", f["href"])][0]
            if "https:" not in download_link:
                download_link="https:"+download_link
            pdb_links.append(download_link)
        except KeyError as e:
            pass      
    return pdb_links

def find_category_proteins(args):
    '''
    Function to find the links to pdb files for proteins in <category_name>
    '''
    categories, name, site_url, protein_dict = args
    category=categories[name]
    #find all proteins listed in <category> and get their names and respective page links
    proteins=category.find_all(href=re.compile('/motm/\d{3}'))[1::2]
    protein_names=[p.text for p in proteins]
    protein_links=[site_url+p["href"] for p in proteins]
    proteins=dict(zip(protein_names, protein_links))
    
    #for each page in <proteins> find pdb file download links
    for protein in proteins.keys():
        if protein in protein_dict.keys():
            protein_dict[protein]["category"].append(name)
        else:
            #get pdb download links for <protein>
            links=find_links(proteins[protein])
            protein_dict[protein]={
                "links":links,
                "category":[name]
            }
    return protein_dict

def get_protein_links(url):
    '''
    Get links for each protein .pdb file
    ----------
    Parameters
    ----------
        url - str
            url of the pdb website
    '''
    site_url="/".join(url.split("/")[:-2]) #https://pdb101.rcsb.org
    pdb_html=get_html(url)

    category_names=[category.text for category in pdb_html.find_all("a", {"class":"no-underline"})]
    category_html=pdb_html.find_all(id=re.compile('subcategory_\d*'))
    categories=dict(zip(category_names, category_html))

    protein_dict={}
    for cat in tqdm(categories.keys(), desc="Preparing Download"):
        args=[categories, cat, site_url, protein_dict]
        protein_dict=find_category_proteins(args)

    return protein_dict

def download(args):
    '''
    Function to download protein
    ----------
    Parameters 
    ---------
        args - tuple, list
            Contains arguments to download pdb file.
    '''
    protein, protein_dict, dst=args
    protein_name=protein.replace(" ", "-")
    links=protein_dict[protein]["links"]
    save_to=os.path.join(dst, protein_name)
    os.makedirs(save_to, exist_ok=True)
    # print("test")
    for link in links:
        # print("test2")
        fn=link.split("/")[-1]
        fp=os.path.join(save_to, fn)
        urlretrieve(link, fp)

def download_proteins(url, dst):
    '''
    Download pdb files listed in <url>
    ----------
    Parameters
    ----------
        url - str
            url of page containing list of proteins with categories
        dst - str
            path to save proteins to
    '''
    protein_dict=get_protein_links(url)
    # print(protein_dict)
    # args=[(protein, protein_dict, dst) for protein in protein_dict.keys()]

    # with Pool(N_PROCESSES) as pool:
    #     tqdm(pool.imap(download, args), total=len(args))
    #     pool.close()
    #     pool.join()

    for protein in tqdm(protein_dict.keys(), desc="Downloading"):
        args=(protein, protein_dict, dst)
        download(args)
        protein_name=protein.replace(" ", "-")
        links=protein_dict[protein]["links"]
        save_to=os.path.join(dst, protein_name)
        os.makedirs(save_to, exist_ok=True)
        
        for link in links:
            fn=link.split("/")[-1]
            fp=os.path.join(save_to, fn)
            urlretrieve(link, fp)
    
    with open(os.path.join(dst, "proteins.json"), "w") as outfile:
        json.dump(protein_dict, outfile, indent=4, sort_keys=True)
    print("Files saved to %s"%dst)

if __name__=="__main__":
    with open(PARAMS_FP, "r") as file:
        params=json.load(file)

    params=params['scraping-params']
    url=params["url"]
    dst=params["dst"]
    download_proteins(url, dst)
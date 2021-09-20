#-------------------------------------------------------------------------------
#Sept 08, 2021
#Downloads the prerequisite references0 for MetaPro installation
#Also extracts them

import requests
import sys
import os
from datetime import datetime as dt
import time
import shutil

def get_file_with_resume(url, filepath, filename):
    dirname = os.path.dirname(filepath)
    remote_filesize = int(requests.get(url, stream=True).headers['Content-length'])
    if(os.path.isdir(dirname)):
        if(os.path.exists(filepath)):
            local_filesize = os.path.getsize(filepath)
            print(dt.today(), "remote size:", remote_filesize, "local size:", local_filesize)
            if(int(remote_filesize) > int(local_filesize)):
                print(dt.today(), "resume download for:", filename)
                #partial download
                header_string = "bytes=" + str(local_filesize) + "-"
                header={"Range": header_string}
                remaining_bytes = int(remote_filesize) - int(local_filesize)

                with requests.get(url, stream=True, headers=header) as r:
                    r.raise_for_status()
                    with open(filepath, "ab") as file_out:
                        for chunk in r.iter_content(chunk_size=8192):
                            file_out.write(chunk)
                            local_filesize = os.path.getsize(filepath)
                            print(dt.today(), local_filesize, "/", remote_filesize, str(round(100*local_filesize/remote_filesize, 2)) + "%              ", end="\r")
                    print(dt.today(), "finished download:", name)        
                    
            elif(int(remote_filesize) == int(local_filesize)):
                print(dt.today(), "skipping:", filename)

            else:
                print(dt.today(), "local filesize inconsistent with current distribution, deleting and redownloading")
                if(os.path.exists(filepath)):
                    os.remove(filepath)
                    print(dt.today(), "new file download:", filename)
                    with requests.get(url, stream=True) as r:
                            r.raise_for_status()
                            with open(filepath, "wb") as file_out:
                                for chunk in r.iter_content(chunk_size=8192):
                                    file_out.write(chunk)
                                    local_filesize = os.path.getsize(filepath)
                                    print(dt.today(), local_filesize, "/", remote_filesize, str(round(100*local_filesize/remote_filesize, 2)) + "%              ", end="\r")
                            print(dt.today(), "finished download:", name)
        else:
            #folder exists. file doesn't
            print(dt.today(), "Folder exists, file doesn't: new file download:", filename)
            with requests.get(url, stream=True) as r:
                    r.raise_for_status()
                    with open(filepath, "wb") as file_out:
                        for chunk in r.iter_content(chunk_size=8192):
                            file_out.write(chunk)
                            local_filesize = os.path.getsize(filepath)
                            print(dt.today(), local_filesize, "/", remote_filesize, str(round(100*local_filesize/remote_filesize, 2)) + "%              ", end="\r")
                    print(dt.today(), "finished download:", name)
    else:   
        #folder doesn't exist. just download. make the dir too
        dirname = os.path.dirname(filepath)
        if(os.path.isfile(dirname)):
            print(dt.today(), "dirname exists as path: deleting file")
            os.remove(dirname)
        os.makedirs(dirname)
        print(dt.today(), "new file download:", filename)
        with requests.get(url, stream=True) as r:
                r.raise_for_status()
                with open(filepath, "wb") as file_out:
                    for chunk in r.iter_content(chunk_size=8192):
                        file_out.write(chunk)
                        local_filesize = os.path.getsize(filepath)
                        print(dt.today(), local_filesize, "/", remote_filesize, str(round(100*local_filesize/remote_filesize, 2)) + "%              ", end="\r")

                print(dt.today(), "finished download:", name)
           
def start_download(url, path, name, things_to_skip, bypass_log_file):
    print(dt.today(), "starting download:", name)
    get_file_with_resume(url, path, name)
    
    
    if(not name in things_to_skip):
        print(dt.today(), "extracting", name)
        extract_path = os.path.dirname(path).split("/")
        real_extract_path = extract_path[0]
        for item in extract_path[1:-1]:
            real_extract_path = os.path.join(real_extract_path, item)
        bypass_log_file.write(name +"\n")
        shutil.unpack_archive(path, outdir)
        print(dt.today(), "finished extracting", name)
    else:
        print(dt.today(), "skipping extraction:", name)
    print("\n")

def import_bypass_log(bypass_log_path):
    things_to_skip = list()
    if(os.path.exists(bypass_log_path)):
        with open(bypass_log_path, "r") as bypass_log:
            for line in bypass_log:
                cleaned_item = line.strip("\n")
                things_to_skip.append(cleaned_item)
    return things_to_skip


if __name__ == "__main__":
    print("-----------------------------------------------------")
    print("MetaPro reference library downloader v1.0.0")
    
    files_to_download = list()
    #files_to_download.append("all")


    outdir = sys.argv[1] # output dir. 
    if(outdir == "."):
        outdir = os.getcwd()
    if(not os.path.isdir(outdir)):
        print(dt.today(), "new directory. creating")
        os.makedirs(outdir)
        
    bypass_log_file = ""
    bypass_log_path = os.path.join(outdir, "bypass_log.txt")
    things_to_skip = import_bypass_log(bypass_log_path)
    if(os.path.exists(bypass_log_path)):
    	bypass_log_file = open(bypass_log_path, "a")
    else:
    	bypass_log_file = open(bypass_log_path, "w")
    	
        
    print(dt.today(), "downloading all MetaPro libraries to:", outdir)
    print("len sys.argv:", len(sys.argv))
    if(len(sys.argv) > 2):
        count = 2
        for item in sys.argv[2:len(sys.argv)]:
            files_to_download.append(item)
            print("item["+str(count) + "]:", item)
            count+=1
    else:
        print(dt.today(), "defaulting to all-download")
        files_to_download.append("all")

    print(dt.today(), "downloading:", files_to_download)


    if(("detect" in files_to_download) or ("all" in files_to_download)):
        #---------------------------------------------------------
        # Detect V2
        name = "detect_v2_Support_libs"    
        path = os.path.join(outdir, "DETECTv2", "detect_v2.tar.gz")
        url = "https://compsysbio.org/metapro_libs/DETECTv2/detect_v2.tar.gz"
        start_download(url, path, name, things_to_skip, bypass_log_file)
        
    if(("ec_pathway" in files_to_download) or ("all" in files_to_download)):    
        #---------------------------------------------------------
        # EC_pathway
        name = "ec_pathway_file"
        path = os.path.join(outdir, "EC_pathway", "EC_pathway.tar.gz")
        url = "https://compsysbio.org/metapro_libs/EC_pathway/EC_pathway.tar.gz"
        start_download(url, path, name, things_to_skip, bypass_log_file)

    if(("priam" in files_to_download) or ("all" in files_to_download)):
        #---------------------------------------------------------
        # PRIAM 
        name = "priam files"
        path = os.path.join(outdir, "PRIAM_db", "priam.tar.gz")
        url = "https://compsysbio.org/metapro_libs/PRIAM_db/priam.tar.gz"
        start_download(url, path, name, things_to_skip, bypass_log_file)

    if(("rfam" in files_to_download) or ("all" in files_to_download)):
        #---------------------------------------------------------
        #Rfam
        name = "Rfam"
        path = os.path.join(outdir, "Rfam", "rfam.tar.gz")
        url = "https://compsysbio.org/metapro_libs/Rfam/rfam.tar.gz"
        start_download(url, path, name, things_to_skip, bypass_log_file)

    if(("wevote" in files_to_download) or ("all" in files_to_download)):
        #----------------------------------------------------------------
        #Wevote
        name = "wevote"
        path = os.path.join(outdir, "WEVOTE_db", "wevote.tar.gz")
        url = "https://compsysbio.org/metapro_libs/WEVOTE_db/wevote.tar.gz"
        start_download(url, path, name, things_to_skip, bypass_log_file)

    if(("accession2taxid" in files_to_download) or ("all" in files_to_download)):
        #-------------------------------------------------------
        #accession2taxid
        name = "accession2taxid"
        path = os.path.join(outdir, "accession2taxid", "accession2taxid.tar.gz")
        url = "https://compsysbio.org/metapro_libs/accession2taxid/accession2taxid.tar.gz"
        start_download(url, path, name, things_to_skip, bypass_log_file)

    if(("centrifuge" in files_to_download) or ("all" in files_to_download)):
        #--------------------------------------------
        #centrifuge
        name = "centrifuge"
        path = os.path.join(outdir, "centrifuge_db", "nt.tar.gz")
        url = "https://compsysbio.org/metapro_libs/centrifuge_db/nt.tar.gz"
        start_download(url, path, name, things_to_skip, bypass_log_file)

    if(("chocophlan_h3" in files_to_download) or ("all" in files_to_download)):
        #-------------------------------------------------------------
        #ChocoPhlan h3
        name = "chocophlan_h3_chunks"
        path = os.path.join(outdir, name, name+".tar.gz")
        url = "https://compsysbio.org/metapro_libs/"+name + "/" + name + ".tar.gz"
        start_download(url, path, name, things_to_skip, bypass_log_file)

        name = "chocophlan_h3_unique"
        path = os.path.join(outdir, name, name+".tar.gz")
        url = "https://compsysbio.org/metapro_libs/"+name + "/" + name + ".tar.gz"
        start_download(url, path, name, things_to_skip, bypass_log_file)

    if(("kaiju" in files_to_download) or ("all" in files_to_download)):
        #-------------------------------------------------------------
        #kaiju
        name = "kaiju_db"
        path = os.path.join(outdir, "kaiju_db", "kaiju.tar.gz")
        url = "https://compsysbio.org/metapro_libs/kaiju_db/kaiju.tar.gz"
        start_download(url, path, name, things_to_skip, bypass_log_file)

    if(("nr" in files_to_download) or ("all" in files_to_download)):
        #-------------------------------------------------------------
        #nr
        name = "nr"
        path = os.path.join(outdir, "nr", "nr.tar.gz")
        url = "https://compsysbio.org/metapro_libs/nr/nr.tar.gz"
        start_download(url, path, name, things_to_skip, bypass_log_file)

    if(("path_to_superpath" in files_to_download) or ("all" in files_to_download)):
        #--------------------------------------------------------------
        #path_to_superpath
        name = "path_to_superpath"
        path = os.path.join(outdir, "path_to_superpath", "pathway_to_superpathway.tar.gz")
        url = "https://compsysbio.org/metapro_libs/path_to_superpath/pathway_to_superpathway.tar.gz"
        start_download(url, path, name, things_to_skip, bypass_log_file)

    if(("swissprot" in files_to_download) or ("all" in files_to_download)):
        #--------------------------------------------------------------
        #swissprot
        name = "swissprot"
        path = os.path.join(outdir, "swiss_prot_db", "swiss_prot_db.tar.gz")
        url = "https://compsysbio.org/metapro_libs/swiss_prot_db/swiss_prot_db.tar.gz"
        start_download(url, path, name, things_to_skip, bypass_log_file)

    if(("trimmomatic" in files_to_download) or ("all" in files_to_download)):
        #--------------------------------------------------------------
        #trimmomatic
        name = "trimmomatic"
        path = os.path.join(outdir, "trimmomatic_adapters", "trimmomatic_adapters.tar.gz")
        url = "https://compsysbio.org/metapro_libs/trimmomatic_adapters/trimmomatic_adapters.tar.gz"
        start_download(url, path, name, things_to_skip, bypass_log_file)

    if(("univec" in files_to_download) or ("all" in files_to_download)):
        #--------------------------------------------------------------------
        #univec
        name = "univec"
        path = os.path.join(outdir, "univec_core", "UniVec_Core.tar.gz")
        url = "https://compsysbio.org/metapro_libs/univec_core/UniVec_Core.tar.gz"
        start_download(url, path, name, things_to_skip, bypass_log_file)
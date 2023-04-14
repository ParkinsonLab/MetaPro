import sys
import zlib
import gzip
import zipfile
import os
import shutil


#will create file if it doesn't exist
z = zipfile.ZipFile("new_file.zip", "a", zipfile.ZIP_DEFLATED)
#z.write(sys.argv[1])
#z.close()

path = sys.argv[1]
print("path:", path)
for root, dirs, files in os.walk(path):
    print("root:", root)
    print("dirs:", dirs)
    print("files:", files)
    print("===============================")
    for file in files:
        z.write(os.path.join(root, file))
    
z.close()
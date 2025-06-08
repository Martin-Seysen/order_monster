import os
import glob
import shutil
import zipfile
from argparse import ArgumentParser



FILES = ["*.txt", "*.g"]
DIRS = ["axis", "baby_axis", "certificates"]
DEL_SUBDIRS = [ "shelve" ]


def path_join(*args):
    while len(args) and not args[0]:
        args = args[1:]
    return os.path.join(*args)



def yield_intermediate_files():
    for dir in DIRS:
        for file in FILES:
             pattern = path_join(dir, file)
             for filename in glob.glob(pattern):
                  yield filename

def yield_intermediate_dirs():
    for dir in DIRS:
        for subdir in DEL_SUBDIRS:
            yield path_join(dir, subdir)



def remove_intermediate_files(verbose = 0):
    for filename in yield_intermediate_files():
        try:
            os.remove(filename)
            if verbose:
                print(f"File {filename} removed")
        except:
            print(f"Could not remove file {filename}")
    for dirname in yield_intermediate_dirs():
        try:
            shutil.rmtree(dirname)
            if verbose:
                print(f"Directory {dirname} removed")
        except:
            if os.path.isdir(dirname):
                print(f"Could not remove directory {dirname}")
            pass


SOURCE_DIRS = ["", "axis", "baby_axis", "utilities", "certificates", "scripts", ]
SOURCE_FILES = ["*.py"]
EXCLUDED = ["zipme.py"]
ADDED_FILES = ["readme.txt"]
ZIP_FILE = "axis_orbits.zip"
PATH_PREFIX = ""



def yield_source_files():
    for dir in SOURCE_DIRS:
        for file in SOURCE_FILES:
             pattern = path_join(dir, file)
             for filename in glob.glob(pattern):
                  ok = True
                  for exclude in EXCLUDED:
                      if filename in glob.glob(exclude):
                          ok = False
                  if ok:
                      yield filename
    for filename in ADDED_FILES:
        yield filename     




def zip_file(zip_filename, prefix, files_to_zip):
    with zipfile.ZipFile(zip_filename, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for file_path in files_to_zip:
            if os.path.isfile(file_path):
                # Preserve the relative path and add the prefix
                arcname = os.path.join(prefix, file_path)
                arcname = arcname.replace("\\", "/")
                zipf.write(file_path, arcname)
                print(f"Added: {file_path} as {arcname}")
            else:
                #print(f"Skipped (not found): {file_path}")
                pass


def main():
    usage = "%prog [options]"
    parser = ArgumentParser(usage)
    parser.add_argument("-r",  dest="remove", action="store_true",
        help="remove all intermediate files and internal tables")
    parser.add_argument("-l",  dest="lst", action="store_true",
        help="list all intermediate files found")
    parser.add_argument("-q",  dest="quiet", action="store_true",
        help="quiet operation" )
    parser.add_argument("-z",  dest="zip", action="store_true",
        help = f"zip source files to the file '{ZIP_FILE}'" )
    
    options  = parser.parse_args()

    if options.lst:
        print("Intermedate files found:")
        for filename in yield_intermediate_files():
            print(" file ", filename)
        for dirname in yield_intermediate_files():
            print(" directory ", dirname)
    if options.remove:
        remove_intermediate_files(not options.quiet)
    if options.zip:
        zip_file(ZIP_FILE, PATH_PREFIX, yield_source_files())


if __name__ == "__main__":
    main()



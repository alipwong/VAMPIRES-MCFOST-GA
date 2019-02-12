import subprocess
import time

def clear_directory(path, sleep = 1):
    '''deletes all the files in a directory'''
    subprocess.Popen("rm -rf {}*".format(path), shell=True)
    time.sleep(sleep)

import subprocess

def clear_directory(path):
    '''deletes all the files in a directory'''
    subprocess.Popen("rm -rf {}*".format(path), shell=True)

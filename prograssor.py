from progressorsClass import *
from time import sleep
pb = progressbar(8, "*")
for count in range(1, 9):
    pb.progress(count)
    sleep(2)
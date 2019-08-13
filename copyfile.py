import os
from shutil import copy2
from glob import glob

for j in range(128):
    os.mkdir('D:/yeongu/plots/8pc_n0_slice/y_%s_slice' % j)
    for i in range(251, 501):
        copy2('D:/yeongu/plots/8pc_n0_slice/%s_slice/n0_slice_%s_%s.png' % (i,i,j),'D:/yeongu/plots/8pc_n0_slice/y_%s_slice' % (j))
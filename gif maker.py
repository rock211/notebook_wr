import imageio
base = 'D/yeongu/plots/proj/'
id = 'ICM1'
dir = 'xproj'

filenames =[]
for i in range(250,500):
    filenames.append(base+id+dir+'/'+'surf_RPS_8pc_ICM1_xproj_%s.png' % i)
print filenames

images = []
for filename in filenames:
    images.append(imageio.imread(filename))
imageio.mimsave('movie.gif', images)



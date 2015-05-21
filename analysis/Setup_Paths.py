import os
try:
    homedir='/Users/jucker/'
    os.listdir(homedir)
except:
    homedir='/Users/mjucker/'

progdir=homedir+'Dropbox/Python/'
repdir =homedir+'Repositories/'
savedir=homedir+'Downloads/'

print 'homedir: ',homedir
print 'repdir:  ',repdir
print 'progdir: ',progdir
print 'savedir: ',savedir

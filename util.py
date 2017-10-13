import os
import glob

def first_last_snap(dirname, stem="snapshot_", fillwidth=4):
	filelist = map(os.path.basename, glob.glob(os.path.join(dirname, stem) + "*"))
	filelist.sort()
	first = int(filelist[0][len(stem):])
	last = int(filelist[-1][len(stem):])

	return first, last
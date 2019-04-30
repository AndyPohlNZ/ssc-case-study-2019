from PIL import Image
import os
import numpy as np
import csv

### read training labels ####
os.chdir("/Users/lq/Downloads/SSC19-master")#please change the directory to your working path
training_image_name = []
training_cell_count = []
training_blur_level = []
training_stain = []
with open('train_label.csv', newline='') as f:
    reader = csv.reader(f)
    next(reader,None)
    for row in reader:
        training_image_name.append(row[1])
        training_cell_count.append(row[2])
        training_blur_level.append(row[3])
        training_stain.append(row[4])

### read training images ####
os.chdir("/Users/lq/Downloads/SSC19-master/train/")#please change the directory to your working path
##create a list (training_list) containing pixel value of 2400 training images
traning_list = []
for i in training_image_name:
	im = Image.open(i)
	#im.show()
	imarray = np.array(im)
	traning_list.append(imarray)



### read test labels ####
os.chdir("/Users/lq/Downloads/SSC19-master")#please change the directory to your working path
testing_image_name = []
testing_cell_count = []
testing_blur_level = []
testing_stain = []
with open('test_label.csv', newline='') as f:
    reader = csv.reader(f)
    next(reader,None)
    for row in reader:
        testing_image_name.append(row[1])
        testing_cell_count.append(row[2])
        testing_blur_level.append(row[3])
        testing_stain.append(row[4])

## read testing images ####        
os.chdir("/Users/lq/Downloads/SSC19-master/test/")#please change the directory to your working path
##create a list (testing_list) containing pixel value of 1200 testing images
testing_list = []
for i in testing_image_name:
	im = Image.open(i)
	#im.show()
	imarray = np.array(im)
	testing_list.append(imarray)











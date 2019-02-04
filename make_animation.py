
import cv2
import os

image_folder = '/home/hag007/Desktop/nn/'
video_name = image_folder+'video.avi'

images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
frame = cv2.imread(os.path.join(image_folder, images[0]))
height, width, layers = frame.shape

video = cv2.VideoWriter(video_name, 0, 5, (width,height))

for image in sorted(images, key=lambda x : int(x.split("_")[-1].split(".")[0])):
    video.write(cv2.imread(os.path.join(image_folder, image)))

cv2.destroyAllWindows()
video.release()
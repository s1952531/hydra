#!/bin/bash

# Loop through numbers from 0.0 to 40.0 with step 0.1, adjust as desired
for i in $(seq 0.0 0.1 40.0)
do
  # Use printf to simulate input "l" and the loop number $i
  printf "\n$i\n\n" | python image2.py
  mv l_t$i.png l_t$(printf "%03d" $(echo "$i*10" | bc | cut -d'.' -f1)).png
done
ffmpeg -framerate 10 -i l_t%03d.png -c:v libx264 -pix_fmt yuv420p video.mp4
rm l_t*.png

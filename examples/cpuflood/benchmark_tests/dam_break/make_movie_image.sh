#!/bin/bash -e

ffmpeg -r 10 -i frame%04dfig1.png -c:v libx264 -pix_fmt yuv420p frame.mp4
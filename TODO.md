June 4, 2025:
Trying to fit an estimated anti-clothoid to an actual anti-clothoid curve.
We are fitting the curve doing least squares on the arc-length / RoC^2 plot.
This sorta works alright, but I didn't plot this. However, the location of the
start and alignment of the curve is not as optimal. We are aligning to the first
point and first discrete angle.

Some of the fits return nan, so I was looking at the PrimitiveFitUtils.cpp file.

I am on pause because of the Face Reparam paper.

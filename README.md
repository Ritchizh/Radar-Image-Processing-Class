# Radar-Image-Processing-Class
<p float="left">
<img src = "https://github.com/Ritchizh/Radar-Image-Processing-Class/blob/master/imgs_for_readme/Figure_1.png" height = 250>
<img src = "https://github.com/Ritchizh/Radar-Image-Processing-Class/blob/master/imgs_for_readme/Figure_2.png" height = 220>
<img src = "https://github.com/Ritchizh/Radar-Image-Processing-Class/blob/master/imgs_for_readme/Figure_3.png" height = 220>
</p>

**d_R.npy** -- a sample data file with microwave holography data (SAR principle)

**methods** -- a class this methods to preprocess radar data and reconstruct a radar image.
The most important method is _def focus(self, d)_ -- it's the backpropagation in Fourier frequency domain from the scanning plane to the object's plane.

For detailed information on radar image reconstruction see the classic work:

_Sheen D.M., McMakin D.L., Hall T.E., “Three-dimensional millimeterwave imaging for concealed weapon detection,” IEEE Trans.
Microwave Theory Tech., vol. 49, no. 9, pp. 1581–1592, Sep. 2001._

or its recent non-destructive testing application in our paper:

_M. Chizh, A. Zhuravlev, V. Razevig and S. Ivashov, "Broadband Microwave Imaging for Foam Insulation Diagnostics," 2018 Progress in Electromagnetics Research Symposium (PIERS-Toyama), Toyama, 2018, pp. 1887-1894. DOI: 10.23919/PIERS.2018.8598093_

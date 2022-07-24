# dmag

Written and tested on Ubuntu 22.04.

Creates a depth map (disparity map) from a stereo pair using optical flow.

This repo is an implementation of

[High Accuracy Optical Flow Estimation Based on a Theory for Warping by Thomas Brox, A. Bruhn, N. Papenberg, J. Weickert](https://lmb.informatik.uni-freiburg.de/Publications/2004/Bro04a/)

and

[Robust Optical Flow Estimation by Javier Sánchez Pérez, Nelson Monzón López, Agustín Salgado de la Nuez](https://www.ipol.im/pub/art/2013/21/).

To create the executable, compile the code in directory "dmag" using "make -f Makefile_g/Makefile_O" and then go into the "main" directory and create the exec using "make".

Test cases are given in the "test" directory.

Info about dmag (theory behind it and how to use it) can be found here:

[Depth Map Automatic Generator (DMAG)](https://3dstereophoto.blogspot.com/2013/04/depth-map-automatic-generator-dmag.html)

[Stereo Matching - Variational Methods](https://3dstereophoto.blogspot.com/2012/06/stereo-matching-variational-methods.html)

[Depth Map Generation using Optical Flow](https://www.dropbox.com/s/mpvbuwcwklfh7wn/optical_flow_dmag.pdf?dl=0)

Dependencies (check the Makefiles):

"common" repo

## fold2Bloch

Unfolding of first-principle electronic band structure obtained with [WIEN2k](http://www.wien2k.at) DFT-(L)APW code

### Contributors:
* Oleg Rubel
* Anton Bokhanchuk
* Elias Assmann
* Sheikh Jamil Ahmed

### New features:
[Jul 2022] Adding a transformation matrix that allows to rotate and scale the supercell relative to the primitive cell (e.g., needed for sqrt(2) cells)

[Feb 2022] Plotting the band structure as a raster image (heat map) similar to ARPES

[May 2020] Full support of spinor wave functions. See this thread https://www.mail-archive.com/wien@zeus.theochem.tuwien.ac.at/msg20053.html for more details. (Special thanks to Peter Blaha for helpful suggestions.)

### Documentation, examples, and support

For installation and tutorials, please refer to [wiki page](https://github.com/rubel75/fold2Bloch/wiki)

Examples of "real life" applications can be found in [Phys. Rev. B **90**, 115202 (2014)](http://olegrubel.mcmaster.ca/publications/2014/Rubel_PRB_90_115202.pdf) and [Phys. Rev. Materials **2**, 114604 (2018)](http://olegrubel.mcmaster.ca/publications/2018/Zheng_PRMat_2_2018.pdf)

Please communicate your feedback, support or feature requests via WIEN2k [mailing list](http://www.wien2k.at/reg_user/mailing_list)

### References
If you find the results useful and publishable, we will appreciate citing the following papers:
* O. Rubel, A. Bokhanchuk, S. J. Ahmed, and E. Assmann "Unfolding the band structure of disordered solids: from bound states to high-mobility Kane fermions", [Phys. Rev. B **90**, 115202 (2014)](http://olegrubel.mcmaster.ca/publications/2014/Rubel_PRB_90_115202.pdf).
* L.-W. Wang, L. Bellaiche, S.-H. Wei, and A. Zunger "Majority representation of alloy electronic states", [Phys. Rev. Lett. **80**, 4725 (1998)](https://doi.org/10.1103/PhysRevLett.80.4725).

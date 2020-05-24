## fold2Bloch

Unfolding of first-principle electronic band structure obtained with [WIEN2k](http://www.wien2k.at) DFT-(L)APW code

### Contributors:
* Anton Bokhanchuk
* Elias Assmann
* Sheikh Jamil Ahmed
* Oleg Rubel

For installation and tutorials, please refer to [wiki page](https://github.com/rubel75/fold2Bloch/wiki)

Examples of "real life" applications can be found in [arXiv:1405.4218](http://arxiv.org/abs/1405.4218)

Please communicate your feedback, support or feature requests via WIEN2k [mailing list](http://www.wien2k.at/reg_user/mailing_list)

Compile: ifort -assume realloc_lhs -assume byterecl -g -traceback -check all -debug all -free util.F line_count.F90 NewK.F90 Sort.F90 SortC.F90 fold2Bloch.F90 -o fold2Bloch

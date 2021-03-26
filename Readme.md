# gr-burst (original readme)

This out of tree module contains a number of blocks which help enable building simple burst PSK modems in GNU Radio

install python-bitarray and build/install gr-mapper

----------------------------------------------------

## Dependencies

https://github.com/osh/gr-eventstream  
https://github.com/gr-vt/gr-mapper  
https://github.com/osh/gr-pyqt
https://github.com/sandialabs/gr-pdu_utils (optional)

## Building
>This module requires **Gnuradio 3.7.x**.
>Build is pretty standard:
```
mkdir build
cd build
cmake ..
make
sudo make install
sudo ldconfig
```
## Building local

>**Gnuradio 3.7.x** is installed to `$HOME/gr3.7`:

```
cd gr3.7
source gr3.7-source.env
cd ..
cd gr-gsSDR
mkdir build 
cd build
cmake ../ -Wno-dev -DCMAKE_INSTALL_PREFIX=~/gr3.7.13.4 
make install
sudo ldconfig
```

## gr-burst extension

>  **QPSK burst synchronization block** has been extended and improved:

- optional **signal decimation after coarse frequency synchronization** (we can compensate bigger freq. offset)
- for example: `4 samples per symbol extends possible freq. compensation window twice in comparison to 2 sps (default)`
- decimation can be enabled/disabled, filter taps can be set ...
- **fixed burst size** (before coarse freq. synchronization - e.g 8192 sample) which means global inicialization of `FFTW objects`
- this option dramatically improves speed of processing !!!
- added **fine frequency synchronization option using second order PLL** (better performance - follow Matlab simulation blocks `qpskSecondOrderPLL.m` and `freq_sync_pll_2_order.m`)
-  **pmt port debug** - enabled / disabled  

![Block structure](https://github.com/pavelfpl/gr-burst/blob/master/qpsk_burst_sync.png)

## Links

https://arxiv.org/pdf/1604.08397.pdf  
https://kirankarra.wordpress.com/machine-learning/qpsk-burst-receiver-synchronization/

### gr-burst and gr-mapper maint-3.8

https://github.com/flynn378/gr-burst  
https://github.com/myersw12/gr-mapper  


# gr-burst (original readme)

This out of tree module contains a number of blocks which help enable building simple burst PSK modems in GNU Radio

install python-bitarray and build/install gr-mapper

----------------------------------------------------

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

### gr-burst maint-3.8

https://github.com/flynn378/gr-burst



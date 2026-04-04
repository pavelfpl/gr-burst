/* -*- c++ -*- */
/* 
 * Copyright 2015/2021 Free Software Foundation, Inc
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifndef INCLUDED_BURST_SYNCHRONIZER_V4_IMPL_H
#define INCLUDED_BURST_SYNCHRONIZER_V4_IMPL_H

#include <gnuradio/fft/fft.h>
#include <burst/synchronizer_v4.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_linalg.h>
#include <gnuradio/filter/fft_filter.h>
#include <mapper/constellation.h>
#include <gnuradio/filter/fir_filter.h>
#include <chrono> 

namespace gr {
  namespace burst {
      
    using namespace std::chrono;   

    class synchronizer_v4_impl : public synchronizer_v4{
     private:
       // Private function
       inline void clearResources();
       void enableDebugMode();
       void set_taps(const std::vector<float> taps);   
       inline void shiftFreq(gr_complex* buf, const int bufLen, const double Fs, const double freq, const double tStart);  
       void qpskBurstCFOCorrect(gr_complex* x_in, const int burstSize, const bool doEstimation);
       void conv_prealloc(const gr_complex* const a, const int aLen, const gr_complex* const b, const int bLen, std::vector<gr_complex> &result);
       void conv_eq_prealloc(const gr_complex* const a, const int aLen, const gr_complex* const b, const int bLen, std::vector<gr_complex> &result);
       void determineOptimalFilter(gsl_vector_complex* w, const gr_complex* const x, const int xLen);
       void toeplitz(const gr_complex* const col,const int M, const gr_complex* const row, const int N, gsl_matrix_complex* T);
       void qpskFirstOrderPLL(const gr_complex* const x, int size, float alpha, gr_complex* y);
       void qpskSecondOrderPLL(const gr_complex* const x, int size, float alpha, float beta, gr_complex* y);
       void handler(pmt::pmt_t msg);

       int32_t m_coarse_timeout;
       time_point<high_resolution_clock> m_start_prev;
       double m_cfoEstimate;
      
       size_t d_pad;
       size_t d_group_delay_offset;
       bool d_even_num_taps;
       filter::kernel::fir_filter_ccf d_fir_ccf;
      
       fft::fft_complex *fftEngine;
       fft::fft_complex *fftEngine1;
       fft::fft_complex *fftEngine2;
       fft::fft_complex *ifftEngineA;
    
       fft::fft_complex *fftEngine11;
       fft::fft_complex *fftEngine22;
       fft::fft_complex *ifftEngineAA;
       
       float d_Fs;
       int d_sps;

       int preFFTEngineFFTSize;
       fft::fft_complex preFFTEngine;

       std::vector<gr_complex> preSyms_fliplr_conj;		
       std::vector<gr_complex> preSyms_xR_fliplr_conj;

       int optimalFilterSize;
       std::vector<gr_complex> wOpt_gr;
       gsl_vector_complex* wOpt;

       int preSymsSize;
       int preSymsRateMatchedSize;

       mapper::constellation d_const;

       bool debugMode;
       bool m_port_debug;
       bool m_decim; 
      
       int m_decimation;
       int m_burst_size;
       int m_burst_size_conv;
       int m_pll_type;
     protected:  
      boost::mutex fp_mutex;
     public:
      synchronizer_v4_impl(double Fs, int sps, std::vector<unsigned char> preamble_bits, std::vector<int> sym_mapping, bool decim, int decimation, int burst_size, int pll_type, bool port_debug, const std::vector<float> taps_,const int coarse_timeout);
      ~synchronizer_v4_impl();
      
      bool stop();

      // Where all the action really happens ...
      // ---------------------------------------
      int work(int noutput_items,
	       gr_vector_const_void_star &input_items,
	       gr_vector_void_star &output_items);

    };

  } // namespace burst
} // namespace gr

#endif /* INCLUDED_BURST_SYNCHRONIZER_V4_IMPL_H */


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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "qa_helpers.h"

#include <gnuradio/fft/fft.h>
#include <gnuradio/io_signature.h>
#include "synchronizer_v4_impl.h"
#include <math.h>
#include <stdio.h>

#include <volk/volk_typedefs.h>
#include <volk/volk.h>

#include <algorithm>    // std::reverse
#include <gnuradio/filter/firdes.h>

#include <chrono> 
#include <iostream> 
#include <vector> 

#define CONST_FIRST_ORDER 0x01
#define CONST_SECOND_ORDER 0x02

#define CONST_DEBUG_MESSAGES
#define CONST_INPUT_VECTOR_EXPORT

using namespace std::chrono; 

namespace gr {
  namespace burst {

// On ARM may produce  volk_32fc_magnitude_32f_* optimized function NAN ...     
// Optimization for ARM may help to increase performace ...
#ifdef __ARM_NEON
    static inline void volk_32fc_magnitude_32f_generic_static(float* magnitudeVector,
                                                   const lv_32fc_t* complexVector,
                                                   unsigned int num_points);  

    static inline void volk_32fc_magnitude_32f_generic_static(float* magnitudeVector,
                                                   const lv_32fc_t* complexVector,
                                                   unsigned int num_points){
        
            const float* complexVectorPtr = (float*)complexVector;
            float* magnitudeVectorPtr = magnitudeVector;
            unsigned int number = 0;
            
            for (number = 0; number < num_points; number++) {
                 const float real = *complexVectorPtr++;
                 const float imag = *complexVectorPtr++;
                 
                 if(real == 0.0 and imag == 0.0){
                    *magnitudeVectorPtr++ = 0.0;
                 }else{
                    *magnitudeVectorPtr++ = sqrtf((real * real) + (imag * imag));
                 }
            }
     }
#endif /* __ARM_NEON */
    
    synchronizer_v4::sptr
    synchronizer_v4::make(double Fs, int sps, std::vector<unsigned char> preamble_bits, std::vector<int> sym_mapping, int decim, int decimation, int burst_size, int pll_type, bool port_debug, const std::vector<float> taps_)
    {
      return gnuradio::get_initial_sptr
        (new synchronizer_v4_impl(Fs, sps, preamble_bits, sym_mapping, decim, decimation, burst_size, pll_type, port_debug, taps_));
    }

    /*
     * The private constructor ...
     * ---------------------------
     */
    synchronizer_v4_impl::synchronizer_v4_impl(double Fs, int sps, std::vector<unsigned char> preamble_bits, std::vector<int> sym_mapping, int decim, int decimation, int burst_size, int pll_type, bool port_debug, const std::vector<float> taps_)
      : gr::sync_block("synchronizer_v4",
              gr::io_signature::make(0, 0, 0),
              gr::io_signature::make(0, 0, 0)),
              optimalFilterSize(48),			                                                 // TODO: the size of the optimal filter should be based on preamble size ... currently hardcoded to 48 taps
              preSymsSize(preamble_bits.size()/2 ),												 // divide by 2 b/c this is QPSK specific, 2 bits/sym
              preSymsRateMatchedSize(preSymsSize*sps),
              d_const(mapper::QPSK, sym_mapping, gr_complex(1,0) ),
              preFFTEngineFFTSize( (int)pow(2,ceil(log2(preamble_bits.size()/2 * sps))) ),		 // 2^nextpow2(preSymsRateMatchedSize)
              preFFTEngine(preFFTEngineFFTSize),												 // initialize the FFT engine to the FFT Size
              d_fir_ccf(1,std::vector<float>(0))
    {
    	message_port_register_in(pmt::mp("cpdus"));
		message_port_register_out(pmt::mp("cpdus"));
		message_port_register_out(pmt::mp("debug_post_cfo"));
		message_port_register_out(pmt::mp("debug_pre_xcorr"));
		set_msg_handler(pmt::mp("cpdus"), boost::bind(&synchronizer_v4_impl::handler, this, _1));
        
        // Set paramaters ...
        // ------------------
    	d_Fs = Fs;
    	d_sps = sps;
       
        m_decim = decim; 
        m_decimation = decimation; if(decim = false) m_decimation = 1;
        m_burst_size = burst_size; if(m_burst_size%2!=0) m_burst_size = m_burst_size-1;
        m_burst_size_conv = (m_burst_size/m_decimation/2)+(2*48)-1;     
        m_port_debug = port_debug;
        m_pll_type = pll_type;
        
        // Default - test - parameters ...
        // -------------------------------
        
        /*
        m_decim = true; 
        m_decimation = 2;
        m_port_debug = false;
        m_burst_size = 8192;  if(m_burst_size%2!=0) m_burst_size = m_burst_size-1;
        m_burst_size_conv = (m_burst_size/m_decimation/2)+(2*48)-1;            
        m_pll_type = 2;
        */
        
        preSyms_fliplr_conj.resize(preSymsSize);
        preSyms_xR_fliplr_conj.resize(preSymsRateMatchedSize);
        wOpt_gr.resize(optimalFilterSize);
    	wOpt = gsl_vector_complex_alloc(optimalFilterSize);

    	// map the preamble bits to symbols, both x1 and x2 ... and flip both vectors for correlation purposes
    	// ---------------------------------------------------------------------------------------------------
        d_const.map(&preamble_bits[0], &preSyms_fliplr_conj[0], preSymsSize, 0);		// after this operation, we are still not flipped and conjugated yet
    	
        // upsample ...
		// ------------
        int jj = 0;
		for(int ii=0; ii<preSymsRateMatchedSize; ii++) {
			if(ii%sps==1) {
				preSyms_xR_fliplr_conj[ii].real(0);
				preSyms_xR_fliplr_conj[ii].imag(0);
			}
			else {
				preSyms_xR_fliplr_conj[ii].real(preSyms_fliplr_conj[jj].real());
				preSyms_xR_fliplr_conj[ii].imag(preSyms_fliplr_conj[jj].imag());
				jj++;
			}
		}
		// TODO: currently the rate matched preamble correlation is not interpolated, just upsampled
		// this seems to work great, but perhaps it is worth doing filtering as well to see if better
		// results can be achieved.  Theoretically, we should get stronger correlation peaks if it is
		// interpolated b/c the intermediary samples would contribute to the correlations strength as well

		std::reverse(preSyms_fliplr_conj.begin(),preSyms_fliplr_conj.end());			// flip the x1 preamble symbols
    	std::reverse(preSyms_xR_fliplr_conj.begin(),preSyms_xR_fliplr_conj.end());		// flip the x2 preamble symbols

    	volk_32fc_conjugate_32fc(&preSyms_fliplr_conj[0], &preSyms_fliplr_conj[0], preSyms_fliplr_conj.size());				// conjugate the x1 preamble symbols
    	volk_32fc_conjugate_32fc(&preSyms_xR_fliplr_conj[0], &preSyms_xR_fliplr_conj[0], preSyms_xR_fliplr_conj.size());	// conjugate the x2 preamble symbols

        std::vector<float> taps;
        
        if(taps_.size()!=0){
           taps = taps_; 
        }else{
           taps = {0.0018798966193571687, -0.0005957852699793875, -0.002913221251219511, -0.002573465695604682, 0.0005119061679579318, 0.003337589092552662, 0.002398500684648752, -0.0023926629219204187, -0.006377795245498419, -0.0037021359894424677, 0.0064182160422205925, 0.016373006626963615, 0.01431208848953247, -0.00553062604740262, -0.03386737033724785, -0.047264810651540756, -0.02122040092945099, 0.05183464661240578, 0.15228667855262756, 0.23982208967208862, 0.27452731132507324, 0.23982208967208862, 0.15228667855262756, 0.05183464661240578, -0.02122040092945099, -0.047264810651540756, -0.03386737033724785, -0.00553062604740262, 0.01431208848953247, 0.016373006626963615, 0.0064182160422205925, -0.0037021359894424677, -0.006377795245498419, -0.0023926629219204187, 0.002398500684648752, 0.003337589092552662, 0.0005119061679579318, -0.002573465695604682, -0.002913221251219511, -0.0005957852699793875, 0.0018798966193571687}; 
        }
        
        // To avoid memory access issues within the filterNdec function of gnuradio ...
        // ----------------------------------------------------------------------------
        // Pad the input at the front and back based on the volk alignemt ...
        d_pad = volk_get_alignment() / sizeof(float);
        set_taps(taps);
        
        fftEngine = new fft::fft_complex(m_burst_size);
        fftEngine1= new fft::fft_complex(m_burst_size_conv);
        fftEngine2 = new fft::fft_complex(m_burst_size_conv);
        ifftEngineA = new fft::fft_complex(m_burst_size_conv,false);
        
        fftEngine11= new fft::fft_complex(2*preSymsSize-1);
        fftEngine22 = new fft::fft_complex(2*preSymsSize-1);
        ifftEngineAA = new fft::fft_complex(2*preSymsSize-1,false);
        
    	debugMode = false;
    }

    /*
     * Our virtual destructor ...
     * --------------------------
     */
    synchronizer_v4_impl::~synchronizer_v4_impl(){
        
       gr::thread::scoped_lock lock(fp_mutex);           // shared resources ...
       
       clearResources();
        
    }
    
    
    // stop - public function ...
    // --------------------------
    bool synchronizer_v4_impl::stop(){
        
        gr::thread::scoped_lock lock(fp_mutex);          // shared resources ...
        
        // Clear resources ...
        clearResources();
        
        return true;
    }
    
    // clearResources - private function ...
    // -------------------------------------
    inline void synchronizer_v4_impl::clearResources(){
        
       if(wOpt!=0) {gsl_vector_complex_free(wOpt); wOpt=0;}
        
       if(fftEngine!=0){ delete fftEngine; fftEngine = 0;}
       
       if(fftEngine1!=0){ delete fftEngine1; fftEngine1 = 0;}
       if(fftEngine2!=0){ delete fftEngine2; fftEngine2 = 0;}
       if(ifftEngineA!=0){ delete ifftEngineA; ifftEngineA = 0;}
       
       if(fftEngine11!=0){ delete fftEngine11; fftEngine11 = 0;}
       if(fftEngine22!=0){ delete fftEngine22; fftEngine22 = 0;}
       if(ifftEngineAA!=0){ delete ifftEngineAA; ifftEngineAA = 0;}
    }
    
    // set_taps - private function ...
    // -------------------------------
    void synchronizer_v4_impl::enableDebugMode() {
    	
        debugMode = true;
    }
    
    // set_taps - private function ...
    // -------------------------------
    void synchronizer_v4_impl::set_taps(std::vector<float> taps){
        
        size_t tap_len = taps.size();
        d_fir_ccf.set_taps(taps);

        if(tap_len % 2) {
            d_group_delay_offset = (tap_len - 1)/2;
            d_even_num_taps = false;
        }else {
            d_group_delay_offset = (tap_len)/2;
            d_even_num_taps = true;
            GR_LOG_WARN(d_logger,
                        "PERFORMANCE IMPACT: Even number of taps requires inefficient manual "
                        "adjustiment of burst time");
        }
    }
    
    // qpskBurstCFOCorrect - private function ...
    // ------------------------------------------
    float synchronizer_v4_impl::qpskBurstCFOCorrect(gr_complex* x_in, int burstSize) {
		
        // TODO/done: maybe we can have a fixed burst size, only calculate cfo over a certain burst size
		// then we can statically create the fft engine in the constructor, will probably be a lot faster
		
        /*
        if(burstSize%2!=0) {
			burstSize = burstSize-1;
		}
		*/
		
		// Using global fftEngine object create in constructor - see TODO/done ... 
		gr_complex* fftInBuf = fftEngine->get_inbuf();
		gr_complex* fftOutBuf = fftEngine->get_outbuf();
		std::vector<float> fftOutAbs(burstSize, 0.0);
		
        // Fill the FFT buffer ...
		volk_32fc_s32f_power_32fc(fftInBuf, x_in, 4, burstSize);		// compute signal^4 FFT (b/c this is QPSK specific)

		if(debugMode) {
			std::string filename = "/tmp/gr_cfoFFTInput.txt";
			std::vector<gr_complex> b(burstSize);
			b.assign(fftInBuf, fftInBuf+burstSize);
			qa_helpers::writeComplexFile(filename, b);
		}

		// Compute the fft ...
		fftEngine->execute();
		
        // Take absolute value of fft output and find argmax simultaneously ...
		int maxIdx = 0;
		float maxVal = -100000;
		for(int ii=0; ii<burstSize; ii++) {
			fftOutAbs[ii] = pow(fftOutBuf[ii].real(), 2) + pow(fftOutBuf[ii].imag(), 2);
			
            if(fftOutAbs[ii]>maxVal) {
				maxVal = fftOutAbs[ii];
				maxIdx = ii;
			}
		}

		if(debugMode) {
			std::string filename = "/tmp/gr_cfoFFTOut_abs.txt";
			qa_helpers::writeFloatFile(filename, fftOutAbs);
		}

		// Account for fftshift ...
		if(maxIdx<burstSize/2) {
			maxIdx = maxIdx + burstSize/2;
		}
		else {
			maxIdx = maxIdx - burstSize/2;
		}
		
		double cfoEstimate = (-d_Fs/2.0+maxIdx*(d_Fs)/(burstSize-1)) / 4.0;

		if(debugMode) {
			std::string filename = "/tmp/gr_cfoEstimate.txt";
			std::vector<float> cfoEstimateVec(1);
			cfoEstimateVec[0] = cfoEstimate;
			qa_helpers::writeFloatFile(filename, cfoEstimateVec);
		}

		// Correct for cfo ...
		// -------------------
		shiftFreq(&x_in[0], burstSize, d_Fs, cfoEstimate, 0);

        return cfoEstimate;
	}

    // qpskBurstCFOCorrect - private (inline) function ...
    // ---------------------------------------------------
    inline void synchronizer_v4_impl::shiftFreq(gr_complex* buf, int bufLen, double Fs, double freq, double tStart) {
		
        std::complex<float> phInc(cos(2*M_PI*-freq*1.0/Fs),sin(2*M_PI*-freq*1.0/Fs));
		std::complex<float> phStart(cos(2*M_PI*-freq*tStart),sin(2*M_PI*-freq*tStart));
		volk_32fc_s32fc_x2_rotator_32fc(buf, buf, phInc, &phStart, bufLen);
    }
    
    // conv - private function ...
    // ---------------------------
     void synchronizer_v4_impl::conv(gr_complex* a, int aLen, const gr_complex* b, int bLen, std::vector<gr_complex> &result){
    	int fftSize = aLen+bLen-1;
		result.resize(fftSize);

		fft::fft_complex fftEngine1 = fft::fft_complex(fftSize);
		gr_complex* fftInBuf = fftEngine1.get_inbuf();
		memcpy(fftInBuf, a, sizeof(gr_complex)*aLen);
		memset(fftInBuf+aLen, 0, sizeof(gr_complex)*(fftSize-aLen));
		fftEngine1.execute();
		gr_complex* fft1OutBuf = fftEngine1.get_outbuf();

		fft::fft_complex fftEngine2 = fft::fft_complex(fftSize);
		fftInBuf = fftEngine2.get_inbuf();
		memcpy(fftInBuf, b, sizeof(gr_complex)*bLen);
		memset(fftInBuf+bLen, 0, sizeof(gr_complex)*(fftSize-bLen));
		fftEngine2.execute();
		gr_complex* fft2OutBuf = fftEngine2.get_outbuf();

		// multiply the fft outputs
		fft::fft_complex ifftEngine = fft::fft_complex(fftSize, false);
		gr_complex* ifftInBuf = ifftEngine.get_inbuf();
		volk_32fc_x2_multiply_32fc(ifftInBuf, fft1OutBuf, fft2OutBuf, fftSize);

		// ifft
		ifftEngine.execute();

		// scale back
		volk_32fc_s32fc_multiply_32fc(&result[0], ifftEngine.get_outbuf(), 1.0/fftSize, fftSize);
	}
    

    // conv_eq_prealloc - private function ...
    // ---------------------------------------
    void synchronizer_v4_impl::conv_eq_prealloc(gr_complex* a, int aLen, const gr_complex* b, int bLen, std::vector<gr_complex> &result){
        
        int fftSize = 2*preSymsSize-1;
		result.resize(fftSize);
        
		// Created statically in constructor: fft::fft_complex fftEngine11 = fft::fft_complex(fftSize);
		gr_complex* fftInBuf_0 = fftEngine11->get_inbuf();
		memcpy(fftInBuf_0, a, sizeof(gr_complex)*aLen);
        memset(fftInBuf_0+aLen,0, sizeof(gr_complex)*(fftSize-aLen));
		fftEngine11->execute();
		gr_complex* fft1OutBuf = fftEngine11->get_outbuf();
        
		// Created statically in constructor: fft::fft_complex fftEngine22 = fft::fft_complex(fftSize);
		gr_complex* fftInBuf_1 = fftEngine22->get_inbuf();
        memcpy(fftInBuf_1, b, sizeof(gr_complex)*bLen);
        memset(fftInBuf_1+bLen, 0, sizeof(gr_complex)*(fftSize-bLen));
		fftEngine22->execute();
		gr_complex* fft2OutBuf = fftEngine22->get_outbuf();

		// Multiply the fft outputs ...
		// Created statically in constructor: fft_complex ifftEngine = fft::fft_complex(fftSize, false);
		gr_complex* ifftInBuf = ifftEngineAA->get_inbuf();
		volk_32fc_x2_multiply_32fc(ifftInBuf, fft1OutBuf, fft2OutBuf, fftSize);

		// ifft
		ifftEngineAA->execute();

		// scale back
		volk_32fc_s32fc_multiply_32fc(&result[0], ifftEngineAA->get_outbuf(), 1.0/fftSize, fftSize);
        
	}
	
    // conv_prealloc - private function ...
    // ------------------------------------
    void synchronizer_v4_impl::conv_prealloc(gr_complex* a, int aLen, const gr_complex* b, int bLen, std::vector<gr_complex> &result){
    	
        int fftSize = m_burst_size_conv;
		result.resize(fftSize);
        
#ifdef CONST_DEBUG_MESSAGES
        /*
        std::cout<<"FFT size conv_preamble 1: "<<aLen<<std::endl;
        std::cout<<"FFT size conv_preamble 2: "<<bLen<<std::endl;
        std::cout<<"FFT size conv_preamble (1+2)-1:"<<fftSize<<std::endl;
        */
#endif
		// Created statically in constructor: fft::fft_complex fftEngine11 = fft::fft_complex(fftSize);
		gr_complex* fftInBuf_0 = fftEngine1->get_inbuf();
		memcpy(fftInBuf_0, a, sizeof(gr_complex)*aLen);
        memset(fftInBuf_0+aLen,0, sizeof(gr_complex)*(fftSize-aLen));
		fftEngine1->execute();
		gr_complex* fft1OutBuf = fftEngine1->get_outbuf();
        
		// Created statically in constructor: fft::fft_complex fftEngine22 = fft::fft_complex(fftSize);
		gr_complex* fftInBuf_1 = fftEngine2->get_inbuf();
		memcpy(fftInBuf_1, b, sizeof(gr_complex)*bLen);
		memset(fftInBuf_1+bLen, 0, sizeof(gr_complex)*(fftSize-bLen));
		fftEngine2->execute();
		gr_complex* fft2OutBuf = fftEngine2->get_outbuf();

		// Multiply the fft outputs ...
		// Created statically in constructor: fft_complex ifftEngine = fft::fft_complex(fftSize, false);
		gr_complex* ifftInBuf = ifftEngineA->get_inbuf();
		volk_32fc_x2_multiply_32fc(ifftInBuf, fft1OutBuf, fft2OutBuf, fftSize);

		// ifft
		ifftEngineA->execute();

		// scale back
		volk_32fc_s32fc_multiply_32fc(&result[0], ifftEngineA->get_outbuf(), 1.0/fftSize, fftSize);
	}
	
    // determineOptimalFilter - private function ...
    // ---------------------------------------------
    void synchronizer_v4_impl::determineOptimalFilter(gsl_vector_complex* w, gr_complex* x, int xLen) {
    	
        gr_complex* fftInBuf = preFFTEngine.get_inbuf();
    	
        // Fill the fft input buffer ...
        memcpy(&fftInBuf[0], &x[0], optimalFilterSize*sizeof(gr_complex));
    	memset(&fftInBuf[optimalFilterSize], 0, (preFFTEngineFFTSize-optimalFilterSize)*sizeof(gr_complex));

    	preFFTEngine.execute();
    	gr_complex* fftOutBuf = preFFTEngine.get_outbuf();

    	if(debugMode) {
           std::string filename = "/tmp/gr_dofFFTOutput.txt";
           std::vector<gr_complex> b(preFFTEngineFFTSize);
           b.assign(fftOutBuf, fftOutBuf+preFFTEngineFFTSize);
           qa_helpers::writeComplexFile(filename, b);
		}

    	// Take the output, and store the absolute value in an array ...
    	gr_complex* ifftInBuf = preFFTEngine.get_inbuf();
        
    	for(int ii=0; ii<preFFTEngineFFTSize; ii++) {
            ifftInBuf[ii].real(pow(std::abs(fftOutBuf[ii]), 2));
            ifftInBuf[ii].imag(0);
    	}

    	if(debugMode) {
           std::string filename = "/tmp/gr_dofIFFTInput.txt";
           std::vector<float> b(preFFTEngineFFTSize);
           for(int ii=0; ii<preFFTEngineFFTSize; ii++) {
				b[ii] = ifftInBuf[ii].real();
           }
           qa_helpers::writeFloatFile(filename, b);
		}

    	preFFTEngine.execute();
    	gr_complex* ifftOutBuf = preFFTEngine.get_outbuf();

    	if(debugMode) {
           std::string filename = "/tmp/gr_dofIFFTOutput.txt";
           std::vector<gr_complex> b(preFFTEngineFFTSize);
           b.assign(ifftOutBuf, ifftOutBuf+preFFTEngineFFTSize);
           qa_helpers::writeComplexFile(filename, b);
		}

    	// Generate the row and col vectors for toeplitz matrix creation and scale as necessary ...
    	// ----------------------------------------------------------------------------------------
		std::vector<gr_complex> row(preSymsSize);
		std::vector<gr_complex> col(preSymsSize);
    	
        for(int ii=0; ii<preSymsSize; ii++) {
    		// Downscale by m and the ifft size ...
    		col[ii] = std::conj(ifftOutBuf[ii])*1.0f/(((float)(optimalFilterSize))*((float)(preFFTEngineFFTSize)));
    		row[ii] = ifftOutBuf[ii]*1.0f/(((float)(optimalFilterSize))*((float)(preFFTEngineFFTSize)));
    	}

    	if(debugMode) {
			std::string filename = "/tmp/gr_dofToeplitzCol.txt";
			qa_helpers::writeComplexFile(filename, col);
			std::string filename2 = "/tmp/gr_dofToeplitzRow.txt";
			qa_helpers::writeComplexFile(filename2, row);
		}

    	gsl_matrix_complex* R = gsl_matrix_complex_alloc(preSymsSize, preSymsSize);
    	toeplitz(&col[0], preSymsSize, &row[0], preSymsSize, R);

    	std::vector<gr_complex> xc(preSymsSize+preSymsSize-1);     // -1
        
        // compute correlation between preSyms and x[1:48]
        // -----------------------------------------------
        
    	// the difference in the cross-correlations between the zero-pad and the non-zeropad are small,
        // not sure if we can get away w/ doing no zeropad?? investigate w/ perofrmacne
        
        // conv(&x[0], preSymsSize, &preSyms_fliplr_conj[0], preSymsSize, xc);
        conv_eq_prealloc(&x[0], preSymsSize, &preSyms_fliplr_conj[0], preSymsSize, xc);
        
    	if(debugMode) {
			std::string filename = "/tmp/gr_dofxc.txt";
			qa_helpers::writeComplexFile(filename, xc);

			// Do this to make sure we load P properly as debugging measure ...
			std::vector<gr_complex> P(optimalFilterSize);
			int jj = xc.size()-1;
			for(int ii=0; ii<optimalFilterSize; ii++) {
				P[ii] = gr_complex(xc[jj].real(), -xc[jj].imag() );
				jj--;
			}
			std::string filename2 = "/tmp/gr_P.txt";
			qa_helpers::writeComplexFile(filename2, P);
		}

    	// Make P vector
        gsl_vector_complex* P = gsl_vector_complex_alloc(optimalFilterSize);
    	gsl_complex cval;
    	
        int jj = xc.size()-1;
    	
        for(int ii=0; ii<optimalFilterSize; ii++) {
    		GSL_SET_COMPLEX(&cval, xc[jj].real(), -xc[jj].imag());
    		gsl_vector_complex_set(P, ii, cval);
    		jj--;
    	}

    	// Solve for R ...
    	int s;
    	gsl_permutation * p = gsl_permutation_alloc(optimalFilterSize);
    
        // Sizes ...
    	// ---------
        /*
        if optimalFilterSize = 48 ... the sizes below are shown as such for readability
        R = [ 48x48 ]
        p = [ 48x1 ]
        s = [ 1x1 ]
        P = [ 48x1 ]
        w (wOpt) = [48x1]
        */
        
    	gsl_linalg_complex_LU_decomp (R, p, &s);
    	gsl_linalg_complex_LU_solve (R, p, P, w);

    	// Rescale w happens when we copy w back to a gr vector ...

    	// Free gsl memory ...
    	gsl_matrix_complex_free(R);
		gsl_vector_complex_free(P);
		gsl_permutation_free(p);
    }

    // toeplitz - private function ...
    // -------------------------------
    void synchronizer_v4_impl::toeplitz(gr_complex* col, int M, gr_complex* row, int N, gsl_matrix_complex* T) {
    	
        // Main diagonal and below the main diagonal ...
        gsl_complex colVal, rowVal;
		for( int d=0; d<M; ++d ) {
			for( int i=0; i<M-d; ++i ) {
				GSL_SET_COMPLEX(&colVal, col[d].real(), col[d].imag());
				gsl_matrix_complex_set(T,i+d,i,colVal);
			}
		}

		// Above the main diagonal ...
		for( int d=1; d<N; ++d ) {
			for( int i=0; i<N-d; ++i ) {
				GSL_SET_COMPLEX(&rowVal, row[d].real(), row[d].imag());
				gsl_matrix_complex_set(T,i,i+d,rowVal);
			}
		}
    }

   // qpskFirstOrderPLL - private function ...
    // ----------------------------------------
    void synchronizer_v4_impl::qpskFirstOrderPLL(gr_complex* x, int size, float alpha, gr_complex* y){
    	
        gr_complex phiHat = gr_complex(1,0);
    	gr_complex xHat, er, phiHatT;
	    for(int ii=0; ii<size; ii++) {
			// correct w/ estimated phase
            y[ii] = x[ii]*phiHat;

            // demodulating circuit
			if(y[ii].real()>=0 && y[ii].imag()>=0) {
				xHat.real(M_SQRT1_2);
				xHat.imag(M_SQRT1_2);
			}
			else if(y[ii].real()>=0 && y[ii].imag()<0) {
				xHat.real(M_SQRT1_2);
				xHat.imag(-M_SQRT1_2);
			}
			else if(y[ii].real()<0 && y[ii].imag()<0) {
				xHat.real(-M_SQRT1_2);
				xHat.imag(-M_SQRT1_2);
			}
			else {
				xHat.real(-M_SQRT1_2);
				xHat.imag(M_SQRT1_2);
			}

			// Loop filter to update phase estimate - !!! div with zero !!!
			// ------------------------------------------------------------
			er = std::conj(xHat)*y[ii];
            
            if(std::abs(er) == 0)
               phiHatT = 0;
            else
               phiHatT = er/std::abs(er);
			phiHat = std::conj( std::pow( phiHatT, alpha)) * phiHat;
	    }
    }
    
    // qpskSecondOrderPLL - private function ...
    // -----------------------------------------
    void synchronizer_v4_impl::qpskSecondOrderPLL(gr_complex* x, int size, float alpha, float beta, gr_complex* y) {
            
            float dI = 0, dQ = 0, q = 0;
            float buff_1 = 0.0, buff_2 = 0.0, buff_3 = 0.0;
            
            // better (fine) frequency correction ...
            for(int ii=0; ii<size; ii++) {
               y[ii] = gr_complex(x[ii].real()*cos(buff_3) - x[ii].imag()*sin(buff_3),
                                  x[ii].real()*sin(buff_3) + x[ii].imag()*cos(buff_3)); 
                
               if(y[ii].real() > 0) dI=1.0; else dI=-1; 
               if(y[ii].real() == 0) dI=0.0;
               
               if(y[ii].imag() > 0) dQ=1.0; else dQ=-1; 
               if(y[ii].imag() == 0) dQ=0.0;
               
               q=y[ii].imag()*dI-y[ii].real()*dQ;
               
               buff_1=buff_1+q*beta;
               buff_2=buff_1+q*alpha;
               buff_3=buff_3-buff_2;  
               
               // -- no carrier - only phase 2*pi*fc/fs --
               if(buff_3>2*M_PI) buff_3=buff_3-2*M_PI;
          }
    }

    // handler - private function ...
    // ------------------------------
    void synchronizer_v4_impl::handler(pmt::pmt_t msg){
        
        gr::thread::scoped_lock lock(fp_mutex);           // shared resources ...
        
        auto start = high_resolution_clock::now(); 
        
		// 1] Get parameters from the pdu ...
        // ----------------------------------
        
        std::vector<gr_complex> eqBurst(pmt::c32vector_elements(pmt::cdr(msg)));
        pmt::pmt_t meta = pmt::car(msg);
        
        if(eqBurst.size() >  m_burst_size){
           return; // Burst size dimension not valid ...   
        }
        
        eqBurst.resize(m_burst_size);
        
#ifdef CONST_INPUT_VECTOR_EXPORT       
        std::string filename = "/tmp/gr_inputData.txt";
        qa_helpers::writeComplexFile(filename, eqBurst);
#endif
        
		// 2 ] Perform cfo correction ...
		// ------------------------------
		float cfo = 0;
		cfo = qpskBurstCFOCorrect(&eqBurst[0], eqBurst.size());
        
#ifdef CONST_DEBUG_MESSAGES
        std::cout<<"Coarse freq. offset estimate: "<<cfo<<std::endl;
#endif
        
		// 3] FIR filter - option - decimate after COARSE freq correction ...
        // ------------------------------------------------------------------
        
        if(m_decim){
            size_t vlen_in = eqBurst.size();                  // Test input vector size ...
        
            if (vlen_in <= d_fir_ccf.ntaps()) return;         // Not enough data to process ...
        
            size_t vlen_out(d_even_num_taps ? vlen_in + 1 : vlen_in);

            std::vector<gr_complex> d_in(vlen_in + 2 * d_pad + 2 * d_group_delay_offset, 0);
            d_in.insert(
                d_in.begin() + d_pad + d_group_delay_offset, eqBurst.begin() , eqBurst.end());
            
            std::vector<gr_complex> eqIn_decimated_0(vlen_out / m_decimation);

            d_fir_ccf.filterNdec(
                eqIn_decimated_0.data(), d_in.data() + d_pad, eqIn_decimated_0.size(), m_decimation); // Do main FIR filtering ...
            
            eqBurst = eqIn_decimated_0; 
            
#ifdef CONST_INPUT_VECTOR_EXPORT   
            std::string filename = "/tmp/gr_inputDataDecimated.txt";
            qa_helpers::writeComplexFile(filename, eqBurst);
#endif

        }

		if(debugMode) {
           std::string filename = "/tmp/gr_burstCFOCorrected.txt";
           qa_helpers::writeComplexFile(filename, eqBurst);
		}

        // A] publish debug port #1 / should be enabled ...
		// ------------------------------------------------
		if(m_port_debug){
           pmt::pmt_t cfo_vec = pmt::init_c32vector(eqBurst.size(), &eqBurst[0]);
		   message_port_pub(pmt::mp("debug_post_cfo"), pmt::cons( meta, cfo_vec));   
        }
		
        // 4] search for start index/timing ...
		// ------------------------------------
		// division with two - searching only in the first half burst samples now ... 
        std::vector<gr_complex> preCrossCorr_cmplx(preSymsRateMatchedSize+eqBurst.size()/2-1); 
    
        // correlation - using 2xFFT/1xIFFT
		conv_prealloc(&eqBurst[0], eqBurst.size()/2, &preSyms_xR_fliplr_conj[0], preSymsRateMatchedSize, preCrossCorr_cmplx);
        
		std::vector<float> preCrossCorr(preCrossCorr_cmplx.size());
		int maxIdx = 0;
		float maxVal = -99999;
		for(int ii=0; ii<preCrossCorr.size(); ii++) {
			preCrossCorr[ii] = std::abs(preCrossCorr_cmplx[ii]);
			if(preCrossCorr[ii]>maxVal) {
				maxIdx = ii;
				maxVal = preCrossCorr[ii];
			}
		}
		
		int preambleIdxStart = maxIdx - preSymsRateMatchedSize + 1;
		
        if(debugMode) {
			std::string f1 = "/tmp/gr_preCrossCorr.txt";
			qa_helpers::writeFloatFile(f1, preCrossCorr);

			std::string filename = "/tmp/gr_preCrossCorrIdx.txt";
			std::vector<float> preambleIdxStartVec(1);
			preambleIdxStartVec[0] = preambleIdxStart;
			qa_helpers::writeFloatFile(filename, preambleIdxStartVec);
		}

        // B] publish debug port #2 ...
		// ----------------------------
		if(m_port_debug){
		   pmt::pmt_t xcorr_vec = pmt::init_f32vector(preCrossCorr.size(), &preCrossCorr[0]);
		   message_port_pub(pmt::mp("debug_pre_xcorr"), pmt::cons( meta, xcorr_vec ) );
        }
        
		if(preambleIdxStart<0) {
			// -- Means we didn't find the preamble, return  --
			return;
		}
		
		// 5] Decimate the signal - zero_padding ...
        // -----------------------------------------
		int eqIn_decimated_len = (eqBurst.size() - preambleIdxStart + 1)/2;
		std::vector<gr_complex> eqIn_decimated(eqIn_decimated_len, gr_complex(0,0));
		
		for(int ii=0; ii<eqIn_decimated.size(); ii++) {
			eqIn_decimated[ii] = eqBurst[preambleIdxStart];
			preambleIdxStart+=2; 
		}
		
        if(debugMode) {
           std::string filename = "/tmp/gr_eqInDecimated.txt";
           qa_helpers::writeComplexFile(filename, eqIn_decimated);
        }

		// 6] Determine optimal filter ...
		// -------------------------------
		determineOptimalFilter(wOpt, &eqIn_decimated[0], eqIn_decimated.size());

		// std::vector<gr_complex> whFilt(eqIn_decimated.size()+optimalFilterSize-1); 
		std::vector<gr_complex> whFilt(m_burst_size_conv);
        
        float maxWopt = -999999;
		float tmp;
		for(int ii=0; ii<optimalFilterSize; ii++) {
			gsl_complex c = gsl_vector_complex_get(wOpt, ii);
			wOpt_gr[ii] = gr_complex(GSL_REAL(c), GSL_IMAG(c));
			tmp = std::abs(wOpt_gr[ii]);
			if(tmp>maxWopt) {
				maxWopt = tmp;
			}
		}

		if(debugMode) {
			std::string filename = "/tmp/gr_wOpt.txt";
			qa_helpers::writeComplexFile(filename, wOpt_gr);
		}
		
		std::complex<float> maxWoptCmplx(1.0/maxWopt,0);
		volk_32fc_s32fc_multiply_32fc(&wOpt_gr[0], &wOpt_gr[0], maxWoptCmplx, wOpt_gr.size());

		conv_prealloc(&eqIn_decimated[0], eqIn_decimated.size(), &wOpt_gr[0], wOpt_gr.size(), whFilt);
        
		if(debugMode) {
			std::string filename2 = "/tmp/gr_wOpt_scaled.txt";
			qa_helpers::writeComplexFile(filename2, wOpt_gr);

			std::string filename3 = "/tmp/gr_whFilt.txt";
			qa_helpers::writeComplexFile(filename3, whFilt);
		}

		// 7] Normalize ...
		// ----------------
		std::vector<float> whFiltAbs(whFilt.size());

#ifdef __ARM_NEON
        volk_32fc_magnitude_32f_generic_static(&whFiltAbs[0], &whFilt[0], whFilt.size())
#else
        volk_32fc_magnitude_32f(&whFiltAbs[0], &whFilt[0], whFilt.size());
#endif        
	    
	    float normFactor = 0;
	    volk_32f_accumulator_s32f(&normFactor, &whFiltAbs[0], whFilt.size());
	    normFactor /= whFilt.size();

	    std::complex<float> normFactorCmplx(1.0/normFactor,0);
	    volk_32fc_s32fc_multiply_32fc(&whFilt[0], &whFilt[0], normFactorCmplx, whFilt.size());

		if(debugMode) {
           std::string filename = "/tmp/gr_whFilt_scaled.txt";
           qa_helpers::writeComplexFile(filename, whFilt);
		}

		// 8] Apply pll ...
		// ----------------
		std::vector<gr_complex> phRecoveredSyms(whFilt.size());
		
        float alpha = 0.0;
        float beta = 0.0;
    
        if(m_pll_type == CONST_FIRST_ORDER){
           alpha = 0.0025;
		   qpskFirstOrderPLL(&whFilt[0], whFilt.size(), alpha, &phRecoveredSyms[0]);  
        }
        
        if(m_pll_type == CONST_SECOND_ORDER){
            
           /*
            * --------------------------
            * -- Sample Matlab script --
            * --------------------------
            symbol_rate = 97.65625e3/2;
            oversample = 1; sig_level = 1;

            fs=symbol_rate*oversample; 
            Ts=1/fs;

            % Loop filter parameters
            Bw=150;              % BW of PLL for carrier synchronization
            zeta=1/sqrt(2);      % 1/sqrt(2)
            kp_c=1;              % Proportional constant
            kp_c=kp_c*sig_level; % Proportional constant - account signal level
            k0_c=1;              % k0

            alfa_p1=2*zeta*((2*Bw*Ts)/(zeta+(1/(4*zeta))));
            beta_p2=((2*Bw*Ts)/(zeta+(1/(4*zeta))))^2;

            alfa=alfa_p1*(1/kp_c)*(1/k0_c);
            beta=beta_p2*(1/kp_c)*(1/k0_c);
            */ 
            
            alpha = 0.0082;     // 0.0055; 0.0082;
            beta = 3.3554e-05;  // 1.4913e-05; 3.3554e-05; 
            qpskSecondOrderPLL(&whFilt[0], whFilt.size(), alpha, beta, &phRecoveredSyms[0]);
        }
        
		if(debugMode) {
           std::string filename = "/tmp/gr_phRecoveredSyms.txt";
           qa_helpers::writeComplexFile(filename, phRecoveredSyms);
		}

		// 9] Put into new pdu and send ...
        // --------------------------------
        int offset = (preSymsRateMatchedSize)/2-1;
		pmt::pmt_t newvec = pmt::init_c32vector(phRecoveredSyms.size()-offset, &phRecoveredSyms[offset]);
        meta = pmt::dict_add(meta, pmt::mp("cfo"), pmt::mp(cfo));
        meta = pmt::dict_add(meta, pmt::mp("sync_delay"), pmt::mp(preambleIdxStart));
		msg = pmt::cons(meta, newvec);
		message_port_pub(pmt::mp("cpdus"), msg);
         
        auto stop = high_resolution_clock::now();
		
        // Get duration - to cast it to proper unit use duration cast method ...
        // ---------------------------------------------------------------------
        auto duration = duration_cast<microseconds>(stop - start); 
        
        // Time taken by function ...
        std::cout << "Time taken by function: "<< duration.count()/1000.0 << " miliseconds" << std::endl;
	}

    int
    synchronizer_v4_impl::work(int noutput_items,
			  gr_vector_const_void_star &input_items,
			  gr_vector_void_star &output_items)
    {
        // Do <+signal processing+> ...
        // ----------------------------
        // Tell runtime system how many output items we produced.
        return noutput_items;
    }

  } /* namespace burst */
} /* namespace gr */


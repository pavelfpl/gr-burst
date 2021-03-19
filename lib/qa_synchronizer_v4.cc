/* -*- c++ -*- */
/* 
 * Copyright 2015 Free Software Foundation, Inc
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


#include <gnuradio/attributes.h>
#include <cppunit/TestAssert.h>
#include "qa_synchronizer_v4.h"
#include "qa_helpers.h"

#include <burst/synchronizer_v4.h>

#include <iostream>
#include <fstream>
#include <string>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>

#include <mapper/constellation.h>

namespace gr {
  namespace burst {

    void
    qa_synchronizer_v4::t1()
    {
    	std::cout << "\nRunning SynchronizerV4 Test 1\n";

    	double Fs = 100e3;
    	int sps = 2;
    	unsigned char preamble_bits_arr[96] = {0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0,
								1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0,
								0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
								1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0,
								1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0,
								0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0,
								0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0,
								1, 1, 0, 0, 0};
		std::vector<unsigned char> preamble_bits(&preamble_bits_arr[0], &preamble_bits_arr[0]+96);
		std::vector<int> sym_mapping(4);
		sym_mapping[0] = 0; sym_mapping[1] = 1; sym_mapping[2] = 3; sym_mapping[3] = 2;
        // Add-on parameters ...
        // ---------------------
        bool decim = true; 
        int decimation = 2;
        int burst_size = 8192;             
        int pll_type = 2;
        bool port_debug = false;
        
        std::vector<float> taps = {0.0018798966193571687, -0.0005957852699793875, -0.002913221251219511, -0.002573465695604682, 0.0005119061679579318, 0.003337589092552662, 0.002398500684648752, -0.0023926629219204187, -0.006377795245498419, -0.0037021359894424677, 0.0064182160422205925, 0.016373006626963615, 0.01431208848953247, -0.00553062604740262, -0.03386737033724785, -0.047264810651540756, -0.02122040092945099, 0.05183464661240578, 0.15228667855262756, 0.23982208967208862, 0.27452731132507324, 0.23982208967208862, 0.15228667855262756, 0.05183464661240578, -0.02122040092945099, -0.047264810651540756, -0.03386737033724785, -0.00553062604740262, 0.01431208848953247, 0.016373006626963615, 0.0064182160422205925, -0.0037021359894424677, -0.006377795245498419, -0.0023926629219204187, 0.002398500684648752, 0.003337589092552662, 0.0005119061679579318, -0.002573465695604682, -0.002913221251219511, -0.0005957852699793875, 0.0018798966193571687}; 
        
		synchronizer_v4::sptr sync_v4 =
			synchronizer_v4::make(Fs, sps, preamble_bits, sym_mapping, decim, decimation, burst_size, pll_type, port_debug,taps);
		sync_v4->enableDebugMode();

    	std::string filename = "/tmp/burst1.txt";
    	std::vector<gr_complex> burst1 = qa_helpers::readComplexFile(filename);

		pmt::pmt_t newvec = pmt::init_c32vector(burst1.size(), &burst1[0]);
		pmt::pmt_t msg = pmt::cons( pmt::PMT_NIL, newvec );

		sync_v4->handler(msg);

		double delta = 0.1;

		// check input data
		std::string f1 = "/tmp/gr_inputData.txt";
		std::vector<gr_complex> inputData = qa_helpers::readComplexFile(f1);
		qa_helpers::areComplexVectorsEqual(burst1, inputData, delta);
		std::cout << "Input Data Passed!\n";

		// check burst cfo
		std::string gr_f2 = "/tmp/gr_cfoFFTInput.txt";
		std::vector<gr_complex> gr_cfoFFTInput = qa_helpers::readComplexFile(gr_f2);
		std::string ml_f2 = "/tmp/ml_cfoFFTInput.txt";
		std::vector<gr_complex> ml_cfoFFTInput = qa_helpers::readComplexFile(ml_f2);
		qa_helpers::areComplexVectorsEqual(ml_cfoFFTInput, gr_cfoFFTInput, delta);
		std::cout << "FFT Input Passed!\n";

		std::string gr_f3 = "/tmp/gr_cfoFFTOut_abs.txt";
		std::vector<float> gr_cfoFFTOut_abs = qa_helpers::readFloatFile(gr_f3);
		std::string ml_f3 = "/tmp/ml_cfoFFTOut_abs.txt";
		std::vector<float> ml_cfoFFTOut_abs = qa_helpers::readFloatFile(ml_f3);
		qa_helpers::areFloatVectorsEqual(ml_cfoFFTOut_abs, gr_cfoFFTOut_abs, delta);
		std::cout << "FFT Output Passed!\n";

		std::string gr_f4 = "/tmp/gr_cfoEstimate.txt";
		std::vector<float> gr_cfoEstimate = qa_helpers::readFloatFile(gr_f4);
		std::string ml_f4 = "/tmp/ml_cfoEstimate.txt";
		std::vector<float> ml_cfoEstimate = qa_helpers::readFloatFile(ml_f4);
		qa_helpers::areFloatVectorsEqual(ml_cfoEstimate, gr_cfoEstimate, 0.1);
		std::cout << "CFO Estimate Passed!\n";

		std::string gr_f5 = "/tmp/gr_burstCFOCorrected.txt";
		std::vector<gr_complex> gr_burstCFOCorrected = qa_helpers::readComplexFile(gr_f5);
		std::string ml_f5 = "/tmp/ml_burstCFOCorrected.txt";
		std::vector<gr_complex> ml_burstCFOCorrected = qa_helpers::readComplexFile(ml_f5);
		qa_helpers::areComplexVectorsEqual(ml_burstCFOCorrected, gr_burstCFOCorrected, delta);
		std::cout << "CFO Correction Passed!\n";

		// check preamble cross correlations
		std::string gr_f6 = "/tmp/gr_preCrossCorrIdx.txt";
		std::vector<float> gr_preCrossCorrIdx = qa_helpers::readFloatFile(gr_f6);
		std::string ml_f6 = "/tmp/ml_preCrossCorrIdx.txt";
		std::vector<float> ml_preCrossCorrIdx = qa_helpers::readFloatFile(ml_f6);
		qa_helpers::areFloatVectorsEqual(ml_preCrossCorrIdx, gr_preCrossCorrIdx, 0.1);
		std::cout << "Preamble Start Idx Passed!\n";

		// check eqIn decimated
		std::string gr_f7 = "/tmp/gr_eqInDecimated.txt";
		std::vector<gr_complex> gr_eqInDecimated = qa_helpers::readComplexFile(gr_f7);
		std::string ml_f7 = "/tmp/ml_eqInDecimated.txt";
		std::vector<gr_complex> ml_eqInDecimated = qa_helpers::readComplexFile(ml_f7);
		qa_helpers::areComplexVectorsEqual(ml_eqInDecimated, gr_eqInDecimated, delta);
		std::cout << "EqIn Decimated Passed!\n";

		// check wiener filter calculation
		std::string gr_f8 = "/tmp/gr_dofFFTOutput.txt";
		std::vector<gr_complex> gr_dofFFTOutput = qa_helpers::readComplexFile(gr_f8);
		std::string ml_f8 = "/tmp/ml_dofFFTOutput.txt";
		std::vector<gr_complex> ml_dofFFTOutput = qa_helpers::readComplexFile(ml_f8);
		qa_helpers::areComplexVectorsEqual(ml_dofFFTOutput, gr_dofFFTOutput, delta);
		std::cout << "DOF FFT Output Passed!\n";

		std::string gr_f9 = "/tmp/gr_dofIFFTInput.txt";
		std::vector<float> gr_dofIFFTInput = qa_helpers::readFloatFile(gr_f9);
		std::string ml_f9 = "/tmp/ml_dofIFFTInput.txt";
		std::vector<float> ml_dofIFFTInput = qa_helpers::readFloatFile(ml_f9);
		qa_helpers::areFloatVectorsEqual(ml_dofIFFTInput, gr_dofIFFTInput, delta);
		std::cout << "DOF IFFT Input Passed!\n";

		std::string gr_f10 = "/tmp/gr_dofIFFTOutput.txt";
		std::vector<gr_complex> gr_dofIFFTOutput = qa_helpers::readComplexFile(gr_f10);
		std::string ml_f10 = "/tmp/ml_dofIFFTOutput.txt";
		std::vector<gr_complex> ml_dofIFFTOutput = qa_helpers::readComplexFile(ml_f10);
		qa_helpers::areComplexVectorsEqual(ml_dofIFFTOutput, gr_dofIFFTOutput, delta);
		std::cout << "DOF IFFT Output Passed!\n";

		std::string gr_f11 = "/tmp/gr_dofToeplitzCol.txt";
		std::vector<gr_complex> gr_dofToeplitzCol = qa_helpers::readComplexFile(gr_f11);
		std::string ml_f11 = "/tmp/ml_dofToeplitzCol.txt";
		std::vector<gr_complex> ml_dofToeplitzCol = qa_helpers::readComplexFile(ml_f11);
		qa_helpers::areComplexVectorsEqual(ml_dofToeplitzCol, gr_dofToeplitzCol, delta);
		std::cout << "DOF ToeplitzMat_COL Passed!\n";

		std::string gr_f12 = "/tmp/gr_dofToeplitzRow.txt";
		std::vector<gr_complex> gr_dofToeplitzRow = qa_helpers::readComplexFile(gr_f12);
		std::string ml_f12 = "/tmp/ml_dofToeplitzRow.txt";
		std::vector<gr_complex> ml_dofToeplitzRow = qa_helpers::readComplexFile(ml_f12);
		qa_helpers::areComplexVectorsEqual(ml_dofToeplitzRow, gr_dofToeplitzRow, delta);
		std::cout << "DOF ToeplitzMat_ROW Passed!\n";

		std::string gr_f13 = "/tmp/gr_P.txt";
		std::vector<gr_complex> gr_P = qa_helpers::readComplexFile(gr_f13);
		std::string ml_f13 = "/tmp/ml_P.txt";
		std::vector<gr_complex> ml_P = qa_helpers::readComplexFile(ml_f13);
		qa_helpers::areComplexVectorsEqual(ml_P, gr_P, delta);
		std::cout << "DOF P Passed!\n";

		std::string gr_f14 = "/tmp/gr_wOpt.txt";
		std::vector<gr_complex> gr_wOpt = qa_helpers::readComplexFile(gr_f14);
		std::string ml_f14 = "/tmp/ml_wOpt.txt";
		std::vector<gr_complex> ml_wOpt = qa_helpers::readComplexFile(ml_f14);
		qa_helpers::areComplexVectorsEqual(ml_wOpt, gr_wOpt, 0.1);
		std::cout << "DOF wOpt Passed!\n";

		std::string gr_f15 = "/tmp/gr_wOpt_scaled.txt";
		std::vector<gr_complex> gr_wOpt_scaled = qa_helpers::readComplexFile(gr_f15);
		std::string ml_f15 = "/tmp/ml_wOpt_scaled.txt";
		std::vector<gr_complex> ml_wOpt_scaled = qa_helpers::readComplexFile(ml_f15);
		qa_helpers::areComplexVectorsEqual(ml_wOpt_scaled, gr_wOpt_scaled, 0.001);
		std::cout << "DOF wOpt_scaled Passed!\n";

		std::string gr_f16 = "/tmp/gr_whFilt.txt";
		std::vector<gr_complex> gr_whFilt_unscaled = qa_helpers::readComplexFile(gr_f16);
		std::string ml_f16 = "/tmp/ml_whFilt.txt";
		std::vector<gr_complex> ml_whFilt_unscaled = qa_helpers::readComplexFile(ml_f16);
		qa_helpers::areComplexVectorsEqual(ml_whFilt_unscaled, gr_whFilt_unscaled, 0.001);
		std::cout << "DOF WHFILT UNSCALED Passed!\n";

		std::string gr_f17 = "/tmp/gr_whFilt_scaled.txt";
		std::vector<gr_complex> gr_whFilt = qa_helpers::readComplexFile(gr_f17);
		std::string ml_f17 = "/tmp/ml_whFilt_scaled.txt";
		std::vector<gr_complex> ml_whFilt = qa_helpers::readComplexFile(ml_f17);
		qa_helpers::areComplexVectorsEqual(ml_whFilt, gr_whFilt, 0.001);
		std::cout << "DOF WHFILT Passed!\n";

		std::string gr_f18 = "/tmp/gr_phRecoveredSyms.txt";
		std::vector<gr_complex> gr_phRecoveredSyms = qa_helpers::readComplexFile(gr_f18);
		std::string ml_f18 = "/tmp/ml_phRecoveredSyms.txt";
		std::vector<gr_complex> ml_phRecoveredSyms = qa_helpers::readComplexFile(ml_f18);
		qa_helpers::areComplexVectorsEqual(gr_phRecoveredSyms, ml_phRecoveredSyms, 0.001);
		std::cout << "DOF PLL Passed!\n";

        std::cout << "SynchronizerV4 Test 1 Passed!!\n";
    }

  } /* namespace burst */
} /* namespace gr */


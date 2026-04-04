#!/usr/bin/env python

from gnuradio import gr;
import pmt, sys, pprint, array, math


class framealign(gr.sync_block):
    def __init__(self, skip=24, blocksize=144*10, blocks = 0, shape_bits = False):
        gr.sync_block.__init__(self,"framealign",[],[])
        self.message_port_register_in(pmt.intern("fpdus"));
        self.message_port_register_out(pmt.intern("fpdus"));
        self.set_msg_handler(pmt.intern("fpdus"), self.handler);
        (self.blocksize, self.skip) = (blocksize, skip);
        (self.blocks, self.shape_bits) = (blocks, shape_bits);
        
        print("Block sizes"+str(self.blocksize))
        
    def handler(self, msg):
        meta = pmt.car(msg);
        data_i = pmt.cdr(msg);
        
        # convert pmt -> int list (of bits) -- chop off front and back
        # data = array.array('B', pmt.u8vector_elements(data_i))
        
        data = array.array('f', pmt.f32vector_elements(data_i))
        data = data[self.skip:-1]
        
        nblocks = 0
        
        # Detect n-blocks automatically --> if 0 ...
        if self.blocks == 0:
            nblocks = math.floor(len(data)/self.blocksize); # For mod. with LDPC codes 
        else:
            nblocks = self.blocks                           # For mod. with Reed Solomon and Convolution codes
        
        # IF nblocks == 0 --> return
        if nblocks == 0:
            return
        
        data = data[0:int(nblocks*self.blocksize)];
        
        # Shape bits in case of the mod. with Reed Solomon and Convolution codes ...
        if self.shape_bits:
            
            data_shape = []
                
            for i in range(len(data)): 
                if data[i] > 0:
                    data_shape.append(1.0)  # Soft bit --> 1.0
                else:
                    data_shape.append(-1.0) # Soft bit --> -1.0
            
            # Publish finally - send output ...
            fout = pmt.init_f32vector(len(data_shape), data_shape);
        else:
            # Publish finally - send output ...
            # fout = pmt.init_u8vector(len(data), data.tolist());
            fout = pmt.init_f32vector(len(data), data.tolist());

        # print "framealign: %d in, %d blocks, %d blocksize, %d out"%(len(pmt.f32vector_elements(data_i)), nblocks, self.blocksize, len(data));

        pdu = pmt.cons(meta, fout);
        self.message_port_pub(pmt.intern("fpdus"), pdu);

    def work(self, input_items, output_items):
        pass




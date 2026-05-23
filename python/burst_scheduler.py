#!/usr/bin/env python
from gnuradio import gr
import pmt
import pdulib
class burst_scheduler(gr.sync_block):
    def __init__(self,barrier=10000,slot=1000):
        gr.sync_block.__init__(
            self,
            name="burst_scheduler",
            in_sig = None,
            out_sig = None)

        self.nproduced_val = 0;
        self.message_port_register_in(pmt.intern("sched_pdu"));
        self.message_port_register_in(pmt.intern("nproduced"));
        self.message_port_register_out(pmt.intern("sched_pdu"));
        self.set_msg_handler(pmt.intern("sched_pdu"), self.sched_pdu);
        self.set_msg_handler(pmt.intern("nproduced"), self.nproduced);
        (self.barrier, self.slot) = (barrier, slot);
        self.sched_barrier = 0;

    def sched_pdu(self, pdu):
        sched_time = (self.nproduced_val + self.barrier);                 # Pick a time in the future - default was 10000 ...
        if(sched_time < self.sched_barrier):
            sched_time = self.sched_barrier;
            print "delaying packet to sched barrier"

        sched_time = sched_time - sched_time%self.slot;                   # Round to nearest slot - default was 1000 ...
        event_length = pmt.length(pmt.cdr(pdu));
        # event_length = pmt.to_long(pmt.dict_ref(pmt.car(pdu), pmt.intern("event_length"), pmt.PMT_NIL));
        self.sched_barrier = sched_time + event_length + self.slot;       # Default was 1000 ...
        print "SCHED_EVENT: time=%d, len=%d  "%(sched_time, event_length);

        pdu = pdulib.pdu_arg_add(pdu, pmt.intern("event_time"), pmt.from_uint64(sched_time));
        self.message_port_pub(pmt.intern("sched_pdu"), pdu);
                
        

    def nproduced(self, produced_pmt):
        nproduced = pmt.to_uint64(produced_pmt);
        self.nproduced_val = nproduced;




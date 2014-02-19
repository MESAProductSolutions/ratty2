#!/usr/bin/env python
'''
You need to have KATCP and CORR installed. Get them from http://pypi.python.org/pypi/katcp and http://casper.berkeley.edu/svn/trunk/projects/packetized_correlator/corr-0.4.0/

Hard-coded for 32bit unsigned numbers.
\nAuthor: Jason Manley, Feb 2011.
'''

import corr,time,numpy,struct,sys,logging,ratty2,cal,conf,iniparse,os,valon_synth,socket


class spec:
    def __init__(self, connect=False, log_handler=None, log_level=logging.INFO, **kwargs):

        if log_handler == None: log_handler=corr.log_handlers.DebugLogHandler(100)
            
        self.lh = log_handler
        self.logger = logging.getLogger('RATTY2')
        self.logger.setLevel(log_level)
        self.logger.addHandler(self.lh)

        if not kwargs.has_key('config_file'): kwargs['config_file']='/etc/ratty2/default'

        self.cal=ratty2.cal.cal(**kwargs)
        self.config=self.cal.config

        self.last_acc_cnt=1 #skip the first accumulation (0), which contains junk.

        if connect:
            self.connect()

    def connect(self):
        self.logger.info('Trying to connect to ROACH %s on port %i...'%(self.config['roach_ip_str'],self.config['katcp_port']))
        self.fpga=corr.katcp_wrapper.FpgaClient(self.config['roach_ip_str'],self.config['katcp_port'],timeout=3,logger=self.logger)
        time.sleep(0.2)
        try:
            self.fpga.ping()
            self.logger.info('KATCP connection to FPGA ok.')
        except:
            self.logger.error('KATCP connection failure. Connection to ROACH failed.')
            raise RuntimeError("Connection to FPGA board failed.")

#        self.rf_frontend=corr.katcp_serial.SerialClient(self.config['roach_ip_str'], timeout=3)
#        time.sleep(0.2)
#        try:
#            self.rf_frontend.ping()
#            self.logger.info('KATCP connection to RF frontend ok.')
#        except:
#            self.logger.error('KATCP connection to RF frontend failed.')
#            raise RuntimeError("Connection to RF frontend box failed.")

    def auto_gain(self,print_progress=False):
        """Try to automatically set the RF attenuators. NOT YET IMPLEMENTED FOR RATTY2"""
        #TODO!!!
        raise RuntimeError('Not yet implemented for RATTY2')
        self.logger.info('Attempting automatic RF gain adjustment...')
        if print_progress: print ('Attempting automatic RF gain adjustment...')
        max_n_tries=10
        n_tries=0
        tolerance=1
        rf_gain=self.config['rf_gain_range'][0]
        self.rf_gain_set(rf_gain)
        time.sleep(0.1)
        self.ctrl_set(mrst='pulse',cnt_rst='pulse',clr_status='pulse',flasher_en=True)
        rf_level=self.adc_amplitudes_get()['adc_dbm']
        if self.status_get()['adc_shutdown'] or self.status_get()['adc_overrange']:
            self.logger.error('Your input levels are too high!')
            raise RuntimeError('Your input levels are too high!')

        while (rf_level < self.config['desired_rf_level']-tolerance or rf_level>self.config['desired_rf_level']+tolerance) and n_tries < max_n_tries:
            rf_level=self.adc_amplitudes_get()['adc_dbm']
            difference = self.config['desired_rf_level'] - rf_level
            rf_gain=self.rf_status_get()[1] + difference
            log_str='Gain was %3.1fdB, resulting in an ADC input level of %5.2fdB. Trying gain of %4.2fdB...'%(self.rf_status_get()[1],rf_level,rf_gain)
            self.logger.info(log_str)
            if print_progress: print log_str
            if self.rf_gain < self.config['rf_gain_range'][0]:
                log_str='Gain at minimum, %4.2fdB.'%self.config['rf_gain_range'][0]
                self.logger.warn(log_str)
                if print_progress: print log_str
                self.rf_gain_set(self.config['rf_gain_range'][0])
                break
            elif rf_gain > self.config['rf_gain_range'][1]:
                log_str='Gain at maximum, %4.2fdB.'%self.config['rf_gain_range'][1]
                self.logger.warn(log_str)
                if print_progress: print log_str
                self.rf_gain_set(self.config['rf_gain_range'][1])
                break
            self.rf_gain_set(rf_gain)
            time.sleep(0.1)
            n_tries += 1
        if n_tries >= max_n_tries: 
            log_str='Auto RF gain adjust failed.'
            self.logger.error(log_str)
            if print_progress: print log_str
        else: 
            log_str='Auto RF gain adjust success.'
            if print_progress: print log_str
            self.logger.info(log_str)

    def auto_fft_shift(self):
        """Try to automatically set the FFT shift schedule"""
        self.ctrl_set(mrst='pulse',cnt_rst='pulse',clr_status='pulse',flasher_en=False)
        stat=self.status_get()
        orig_fft_shift = self.config['fft_shift']
        fft_shift_adj = self.config['fft_shift']
        while not(stat['fft_overrange']): 
            fft_shift_adj = fft_shift_adj << 1
            self.fft_shift_set(fft_shift_adj)
            self.ctrl_set(mrst='pulse',cnt_rst='pulse',clr_status='pulse',flasher_en=False)
            stat=self.status_get()
        fft_shift_adj = fft_shift_adj >> 1
        self.fft_shift_set(fft_shift_adj)
        self.config['fft_shift'] = fft_shift_adj & self.config['fft_shift']
        return (fft_shift_adj & self.config['fft_shift'])


#    def rf_band_get(self):
#        """Grabs the current RF switch state. Returns band integer (switch selection)."""
#
#        if (self.rf_frontend.getd(13)==1 and
#            self.rf_frontend.setd(12)==1 and
#            self.rf_frontend.setd(19)==1 and
#            self.rf_frontend.setd(18)==0 and
#            self.rf_frontend.setd(17)==0 and
#            self.rf_frontend.setd(16)==1 and
#            self.rf_frontend.setd(15)==1 and
#            self.rf_frontend.setd(14)==1):
#            return 1
#
#        elif (self.rf_frontend.getd(13)==1 and
#            self.rf_frontend.setd(12)==0 and
#            self.rf_frontend.setd(19)==1 and
#            self.rf_frontend.setd(18)==1 and
#            self.rf_frontend.setd(17)==1 and
#            self.rf_frontend.setd(16)==1 and
#            self.rf_frontend.setd(15)==0 and
#            self.rf_frontend.setd(14)==1):
#            return 2
#
#        elif (self.rf_frontend.getd(13)==1 and
#            self.rf_frontend.setd(12)==1 and
#            self.rf_frontend.setd(19)==0 and
#            self.rf_frontend.setd(18)==1 and
#            self.rf_frontend.setd(17)==1 and
#            self.rf_frontend.setd(16)==0 and
#            self.rf_frontend.setd(15)==1 and
#            self.rf_frontend.setd(14)==1):
#            return 3
#
#        elif (self.rf_frontend.getd(13)==0 and
#            self.rf_frontend.setd(12)==1 and
#            self.rf_frontend.setd(19)==1 and
#            self.rf_frontend.setd(18)==1 and
#            self.rf_frontend.setd(17)==1 and
#            self.rf_frontend.setd(16)==1 and
#            self.rf_frontend.setd(15)==1 and
#            self.rf_frontend.setd(14)==0):
#            return 4
#        else:
#            raise RuntimeError('Invalid RF switch state. RF frontend fault!')
#
    #def rf_band_set(self,band_sel):
    #    """Configures RF switches to select an RF band. Select between 1 and 4.\n
    #    1: 0-828 MHz   2: NA  3: 900-1670 MHz  4: allpass
    #    """
        # RF lineup:
        # 1. RF.IN -- Amp1 --  RF.SW1.4 -- 0-828MHz filterchain -- RF.SW2.3 -- var.att1 -- Amp2 -- var.att2 -- Amp3 -- var.att3 -- Amp4 -- RF.OUT
        # 2.               --  RF.SW1.2 -- N/C                  -- RF.SW2.5 -- 
        # 3.               --  RF.SW1.5 -- 900MHz-1.67GHz chain -- RF.SW2.2 -- 
        # 4.               --  RF.SW1.3 -- N/C                  -- RF.SW2.4 -- 

        #    'switch1_C5' (select port1-port2) : 12,
        #    'switch1_C3' (select port1-port3) : 13,
        #    'switch1_C4' (select port1-port4) : 18,
        #    'switch1_C6' (select port1-port5) : 19,

        #    'switch2_C5' (select port1-port2) : 16,
        #    'switch2_C3' (select port1-port3) : 17,
        #    'switch2_C4' (select port1-port4) : 14,
        #    'switch2_C6' (select port1-port5) : 15,


#        assert band_sel in range(1,5), "Invalid frequence range %i. Valid frequency ranges are %s."%(band_sel,range(1,5))
        # 0 - 828 MHz
#        if band_sel  == 1:
            #self.rf_frontend.setd(13, 1)
            #self.rf_frontend.setd(12, 1)
            #self.rf_frontend.setd(19, 1)
            #self.rf_frontend.setd(18, 0)
            #self.rf_frontend.setd(17, 0)
            #self.rf_frontend.setd(16, 1)
            #self.rf_frontend.setd(15, 1)
            #self.rf_frontend.setd(14, 1)
        # 800 - 1100 MHz - Not implemented
#        elif band_sel == 2:
            #self.rf_frontend.setd(13, 1)
            #self.rf_frontend.setd(12, 0)
            #self.rf_frontend.setd(19, 1)
            #self.rf_frontend.setd(18, 1)
            #self.rf_frontend.setd(17, 1)
            #self.rf_frontend.setd(16, 1)
            #self.rf_frontend.setd(15, 0)
            #self.rf_frontend.setd(14, 1)
        # 900 - 1670 MHz
#        elif band_sel == 3:
            #self.rf_frontend.setd(13, 1)
            #self.rf_frontend.setd(12, 1)
            #self.rf_frontend.setd(19, 0)
            #self.rf_frontend.setd(18, 1)
            #self.rf_frontend.setd(17, 1)
            #self.rf_frontend.setd(16, 0)
            #self.rf_frontend.setd(15, 1)
            #self.rf_frontend.setd(14, 1)
        # currently not connected
#        elif band_sel == 4:
            #self.rf_frontend.setd(13, 0)
            #self.rf_frontend.setd(12, 1)
            #self.rf_frontend.setd(19, 1)
            #self.rf_frontend.setd(18, 1)
            #self.rf_frontend.setd(17, 1)
            #self.rf_frontend.setd(16, 1)
            #self.rf_frontend.setd(15, 1)
            #self.rf_frontend.setd(14, 0)

    def _rf_band_switch_calc(self,rf_band=None):
        """Calculates the bitmap for the RF switches to select an RF band. Select a band between 1 and 4.\n
        0: 0-828 MHz   1: NA  2: 900-1670 MHz  3: NA
        """
        # RF lineup:
        # 1. RF.IN -- Amp1 --  RF.SW1.4 -- 0-828MHz filterchain -- RF.SW2.3 -- var.att1 -- Amp2 -- var.att2 -- Amp3 -- var.     att3 -- Amp4 -- RF.OUT
        # 2.               --  RF.SW1.2 -- N/C                  -- RF.SW2.5 --
        # 3.               --  RF.SW1.5 -- 900MHz-1.67GHz chain -- RF.SW2.2 --
        # 4.               --  RF.SW1.3 -- N/C                  -- RF.SW2.4 --
        #                                       bitOLD: bitNEW:
        #    'switch1_C5' (select port1-port2) : 31,    0
        #    'switch1_C3' (select port1-port3) : 30,    1
        #    'switch1_C4' (select port1-port4) : 29,    2   
        #    'switch1_C6' (select port1-port5) : 28,    3

        #    'switch2_C5' (select port1-port2) : 27,    4
        #    'switch2_C3' (select port1-port3) : 26,    5
        #    'switch2_C4' (select port1-port4) : 25,    6
        #    'switch2_C6' (select port1-port5) : 24,    7
        #rf_bands=[(29,26),(31,24),(28,27),(30,25)]
        rf_bands=[(2,5),(1,7),(3,4),(0,6)]
        if rf_band==None: rf_band=self.config['band_sel']
        assert rf_band<4,"Requested RF band is out of range (0-3)"
        bitmap=(1<<(rf_bands[rf_band][0]))+(1<<(rf_bands[rf_band][1]))
        return rf_band,bitmap

    def fe_set(self,rf_band=None,gain=None):
        """Configures the analogue box: selects bandpass filters and adjusts RF attenuators; updates global config with changes.\n
        Select a band between 0 and 3: \n
        \t 0: 0-828 MHz \n
        \t 1: 750-1100 MHz \n
        \t 2: 900-1670 MHz \n
        \t 3: passthrough \n
        Valid gain range is -94.5 to 0dB; in this case we distribute the gain evenly across the 3 attenuators. \n
        Alternatively, pass a tuple or list to specify the three values explicitly. \n
        If no gain is specified, default to whatever's in the config file \n"""
        rf_band,bitmap=self._rf_band_switch_calc(rf_band=rf_band) #8bits
        self.config['band_sel']=rf_band
        self.logger.info("Selected RF band %i."%rf_band)

        self.config['rf_attens']=self._rf_atten_calc(gain)
        self.config['rf_gain']=self.config['fe_amp']
        self.config['rf_atten']=0.0
        for (att,atten) in enumerate(self.config['rf_attens']):
            bitmap+=(~(int(-atten*2)))<<(8+(6*int(att))) #6 bits each, following on from above rf_band_select.
            self.logger.info("Set attenuator %i to %3.1f"%(att,atten))
            self.config['rf_atten']+=atten
            if self.config['rf_atten_gain_calfiles'][att] != 'none':
                self.config['rf_gain']+=self.get_interpolated_attens(fileName=self.config['rf_atten_gain_calfiles'][att],setpoint=atten)
            else:
                self.config['rf_gain']+=atten
        #print '0x%08X\n'%bitmap
        self.fpga.write_int('rf_ctrl0',bitmap)

    def set_valon(self):        
        HOST = '192.168.14.76'  # The remote host
        PORT = 7148             # The same port as used by the server 
        
        socket._socketobject.read = socket._socketobject.recv
        socket._socketobject.write = socket._socketobject.send
        
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((HOST, PORT))
        
        valon = valon_synth.Synthesizer(s)                
        valon.set_frequency(valon_synth.SYNTH_A,self.config['sample_clk']/1e6)
	time.sleep(3)
	valon.set_frequency(valon_synth.SYNTH_B,self.config['sample_clk']/1e6)
        
    def initialise(self,skip_program=False, clk_check=False, input_sel='Q',print_progress=False):
        """Initialises the system to defaults."""
        if print_progress:
            print '\tProgramming FPGA...',
            sys.stdout.flush()
        if not skip_program:
            self.fpga.upload_program_bof('/etc/ratty2/boffiles/'+self.config['bitstream'],3333)
            #time.sleep(3)
            self.fpga.progdev(self.config['bitstream'])
            if print_progress: print 'ok'
        elif print_progress: print 'skipped'
        
        """Set the Valon synthesizer frequency to sample clock."""
        if print_progress: 
            print '\tConfiguring Valon frequency to %5.1f MHz sample clock...'%(self.config['sample_clk']/1e6),
            sys.stdout.flush()
            
        self.set_valon()
        if print_progress: print 'ok'   
            
        if clk_check: 
            if print_progress:
                print '\tChecking clocks...',
                sys.stdout.flush()
            est_rate=self.clk_check()
            if print_progress: print 'ok, %i MHz'%est_rate
        
        if print_progress:
            print '\tSelecting RF band %i (%i-%i MHz) and adjusting attenuators...'%(self.config['band_sel'],self.config['ignore_low_freq']/1.e6,self.config['ignore_high_freq']/1.e6),
        self.fe_set()
        if print_progress: print 'ok: %3.1fdb total (%2.1fdb, %2.1fdB, %2.1fdB)'%(self.config['rf_atten'],self.config['rf_attens'][0],self.config['rf_attens'][1],self.config['rf_attens'][2])
    
        if print_progress:
            print '\tConfiguring FFT shift schedule...',
            sys.stdout.flush()
        self.fft_shift_set()
        if print_progress: print 'ok'

        if print_progress:
            print '\tConfiguring accumulation period to %4.2f seconds...'%self.config['acc_period'],
            sys.stdout.flush()
        self.acc_time_set(self.config['acc_period'])
        if print_progress: print 'ok'

        if print_progress:
            print '\tClearing status...',
            sys.stdout.flush()
        self.ctrl_set(mrst='pulse',cnt_rst='pulse',clr_status='pulse',flasher_en=False)
        if print_progress: print 'ok'

        stat=self.status_get()
        if stat['adc_shutdown']: 
            log_msg='ADC selfprotect due to overrange!'
            self.logger.error(log_msg)
            if print_progress: print log_msg
        elif stat['adc_overrange']: 
            log_msg='ADC is clipping!'
            self.logger.warn(log_msg)
            if print_progress: print log_msg
        elif stat['fft_overrange']: 
            log_msg='FFT is overflowing!'
            self.logger.error(log_msg)
            if print_progress: print log_msg

    def clk_check(self):
        """Performs a clock check and returns an estimate of the FPGA's clock frequency."""
        est_rate=round(self.fpga.est_brd_clk())
        if est_rate>(self.config['fpga_clk']/1e6 +1) or est_rate<(self.config['fpga_clk']/1e6 -1):
            self.logger.error('FPGA clock rate is %i MHz where we expect it to be %i MHz.'%(est_rate,self.config['fpga_clk']/1e6))
            #raise RuntimeError('FPGA clock rate is %i MHz where we expect it to be %i MHz.'%(est_rate,self.config['fpga_clk']/1e6))
        return est_rate

    def cal_gains(self,low,high):
        """Used as part of system calibration, to derive RF attenuator calibration files. NOT YET IMPLEMENTED FOR RATTY2!"""
        #TODO Implement for RATTY2
        raise RuntimeError('Not yet implmeneted for RATTY2')
        base_gain=self.rf_atten_set((low+high)/2.)
        time.sleep(0.2)
        base_power=self.adc_amplitudes_get()['adc_dbm']-base_gain
        gain_cal=[]
        for g in numpy.arange(low,high+self.config['rf_gain_range'][2],self.config['rf_gain_range'][2]):
            self.rf_atten_set(g)
            time.sleep(0.2)
            gain_cal.append(self.adc_amplitudes_get()['adc_dbm']-base_power)
        return gain_cal

    def get_spectrum(self,last_acc_cnt=None):
        """Gets data from ROACH board and returns the spectra and the state of the roach at the last timestamp.
            Units of 'cal_spectrum' are dBm unless an antenna was specified in your config file, in which case units are dBuV/m.\n
            Performs bandpass correction, fft_scaling adjustment, backs out number of accumulations, RF frontend gain etc.\n"""
        if last_acc_cnt==None: last_acc_cnt=self.last_acc_cnt
        while self.fpga.read_uint('acc_cnt') <= (last_acc_cnt):  #Wait until the next accumulation has been performed. Polling; too bad!
            time.sleep(0.1)
            #print "cnt = " + str(self.fpga.read_uint('acc_cnt'))
        spectrum = numpy.zeros(self.config['n_chans']) 
        for i in range(self.config['n_par_streams']):
            spectrum[i::self.config['n_par_streams']] = numpy.fromstring(self.fpga.read('%s%i'%(self.config['spectrum_bram_out_prefix'],i),self.config['n_chans']/self.config['n_par_streams']*8),dtype=numpy.uint64).byteswap()
        self.last_acc_cnt = self.fpga.read_uint('acc_cnt')
        if self.config['flip_spectrum']: spectrum=spectrum[::-1] 
        stat = self.status_get()
        ampls = self.adc_amplitudes_get()
        if stat['adc_shutdown']: self.logger.error('ADC selfprotect due to overrange!')
        elif stat['adc_overrange']: self.logger.warning('ADC is clipping!')
        elif stat['fft_overrange']: self.logger.error('FFT is overflowing!')
        #print '[%i] %s: input level: %5.2f dBm (ADC %5.2f dBm).'%(last_acc_cnt,time.ctime(timestamp),stat['input_level'],stat['adc_level']),
        #print '\t\tMean raw: %i'%numpy.mean(spectrum)

        cal_spectrum = spectrum
        cal_spectrum /= float(self.config['n_accs'])
        cal_spectrum *= self.config['fft_scale']
        cal_spectrum *= self.config['adc_v_scale_factor']
        #cal_spectrum /= self.config['chan_width']
        #print '\t\tMean adc_v scale: %f'%numpy.mean(cal_spectrum)
        cal_spectrum  = 10*numpy.log10(cal_spectrum)
        #print '\t\tMean log: %f'%numpy.mean(cal_spectrum)
        cal_spectrum -= self.config['rf_gain']
        cal_spectrum -= self.config['pfb_scale_factor']
        cal_spectrum -= self.config['system_bandpass']
        if self.config['antenna_bandpass_calfile'] != 'none':
            cal_spectrum = ratty2.cal.dbm_to_dbuv(cal_spectrum)
            cal_spectrum += self.config['ant_factor']
        #print '\t\tMean ant_bp: %f'%numpy.mean(cal_spectrum)
    

        return {'raw_spectrum':spectrum, 
                'calibrated_spectrum':cal_spectrum,
                'timestamp':time.time(),
                'acc_cnt':self.last_acc_cnt,
                'adc_overrange':stat['adc_overrange'], 
                'fft_overrange': stat['fft_overrange'], 
                'adc_shutdown':stat['adc_shutdown'], 
                'adc_level':ampls['adc_dbm'], 
                'input_level':ampls['input_dbm'], 
                'adc_temp':self.adc_temp_get(), 
                'ambient_temp': self.ambient_temp_get()} 

    def fft_shift_set(self,fft_shift_schedule=None):
        """Sets the FFT shift schedule (divide-by-two) on each FFT stage. 
            Input is an integer representing a binary bitmask for shifting.
            If not specified as a parameter to this function, set the default level from the config file."""
        if fft_shift_schedule==None: fft_shift_schedule=self.config['fft_shift']
        self.fpga.write_int('fft_shift',fft_shift_schedule) 
        self.config['fft_shift']=fft_shift_schedule
        self.config['fft_scale']=2**(cal.bitcnt(fft_shift_schedule))
        self.logger.info("Set FFT shift to %8x (scaling down by %i)."%(fft_shift_schedule,self.config['fft_scale']))

    def fft_shift_get(self):
        """Fetches the current FFT shifting schedule from the hardware."""
        self.config['fft_shift']=self.fpga.read_uint('fft_shift')
        self.config['fft_scale']=2**(cal.bitcnt(self.config['fft_shift']))
        return self.config['fft_shift'] 

    def ctrl_get(self):
        """Reads and decodes the values from the control register."""
        value = self.fpga.read_uint('control')
        return {'mrst':bool(value&(1<<0)),
                'cnt_rst':bool(value&(1<<1)),
                'clr_status':bool(value&(1<<3)),
                'adc_protect_disable':bool(value&(1<<13)),
                'flasher_en':bool(value&(1<<12)),
                'raw':value,
                }

    def ctrl_set(self,**kwargs):
         """Sets bits of all the Fengine control registers. Keeps any previous state.
             \nPossible boolean kwargs:
             \n\t adc_protect_disable 
             \n\t flasher_en
             \n\t clr_status
             \n\t mrst
             \n\t cnt_rst"""

         key_bit_lookup={
             'adc_protect_disable':   13,
             'flasher_en':   12,
             'clr_status':   3,
             'cnt_rst':      1,
             'mrst':         0,
             }
         value = self.ctrl_get()['raw']
         run_cnt=0
         run_cnt_target=1
         while run_cnt < run_cnt_target:
             for key in kwargs:
                 if (kwargs[key] == 'toggle') and (run_cnt==0):
                     value = value ^ (1<<(key_bit_lookup[key]))
                 elif (kwargs[key] == 'pulse'):
                     run_cnt_target = 3
                     if run_cnt == 0: value = value & ~(1<<(key_bit_lookup[key]))
                     elif run_cnt == 1: value = value | (1<<(key_bit_lookup[key]))
                     elif run_cnt == 2: value = value & ~(1<<(key_bit_lookup[key]))
                 elif kwargs[key] == True:
                     value = value | (1<<(key_bit_lookup[key]))
                 elif kwargs[key] == False:
                     value = value & ~(1<<(key_bit_lookup[key]))
                 else:
                     raise RuntimeError("Sorry, you must specify True, False, 'toggle' or 'pulse' for %s."%key)
             self.fpga.write_int('control', value)
             run_cnt = run_cnt +1

    def _rf_atten_calc(self,gain=None,print_progress=False):
        """Calculates the attenuations for each of the 3 RF attenuators in the RF box. \n
        \t Valid range gain is -94.5 to 0dB; in this case we distribute the gain evenly across the 3 attenuators. \n
        \t Alternatively, pass a tuple or list to specify the three values explicitly. \n
        \t If no gain is specified, default to whatever's in the config file \n"""

        #\t Specify 'auto' to attempt automatic gain calibration. \n
        
#        if len is 3, cont
#        elif len is 1, calc 3 attens
#        elif none, pull from config file
#       round to 0.5dB, return.
        

        if type(gain)==list or type(gain)==numpy.ndarray or type(gain)==tuple:
            rf_attens=gain
        elif type(gain)==int  or type(gain)==float:
            rf_attens=[gain/3. for att in range(3)]
        elif gain==None:
            self.logger.info('Using attenuator settings from config file.')
            if self.config.has_key('rf_atten'):
                rf_attens=[self.config['rf_atten']/3. for att in range(3)]
            elif self.config.has_key('rf_attens'):   
                self.logger.info('Using 3 user-supplied attenuator settings from config file.')
            else:
                raise RuntimeError('Unable to figure out your config file\'s frontend gains; please specify rf_atten (float) or rf_attens (list of 3 floats)!')
        else:
            raise RuntimeError('Unable to figure out your requested attenuation; please specify a single float or a list of 3 floats!')
            
        assert len(rf_attens)==3,'Incorect number of gains specified. Please input a list/tuple of 3 numbers'
        for att in range(3):
            assert (-31.5<=rf_attens[att]<=0),"Range for attenuator %i (%3.1f) is out of range (-31.5 to 0)."%(att,rf_attens[att])
        return [round(att*2)/2 for att in rf_attens] 


#        #start by grabbing the defaults from the config file; we'll override these in a minute if the user's asked for something else...
#        #Try to figure out the frontend gain.
#        self.config['rf_gain']=self.config['fe_amp']
#        if self.config.has_key('rf_atten'):
#            self.config['rf_attens']=[]
#            for att in range(3):
#                #TODO: add smarts here.
#                self.config['rf_attens'].append(self.config['rf_atten']/3.)
#
#        if self.config.has_key('rf_attens'):
#            for att in range(3):
#                self.config['rf_attens'][att]=round(self.config['rf_attens'][att]*2)/2
#                if self.config['rf_atten_gain_calfiles'][att] != 'none':
#                    self.config['rf_gain']+=self.get_interpolated_attens(fileName=self.config['rf_atten_gain_calfiles'][att],setpoint=self.config['rf_attens'][att])
#                else:
#                    self.config['rf_gain']+=self.config['rf_attens'][att]
#        else:
#            raise RuntimeError('Unable to figure out your frontend gains; please specify rf_atten (float) or rf_attens (list of 3 floats)!')
#
#
#        self.config['rf_gain']=self.config['fe_amp']
#        set_gain=0
#        if gain==None:
#        #Try to figure out the desired frontend gain from the config file.
#            if self.config.has_key('rf_atten'):
#                self.config['rf_attens']=[]
#                for att in range(3):
#                    #TODO: add smarts here.
#                    self.config['rf_attens'].append(self.config['rf_atten']/3.)
#
#        elif type(gain)==int  or type(gain)==float:
#            for att in range(3):
#                set_gain=-self._rf_atten_write(att,0-gain/3.)   
#                self.config['rf_attens'][att]=set_gain
#
#        elif type(gain)==list or type(gain)==numpy.ndarray or type(a)==tuple:
#            assert len(gain)==3,'Incorect number of gains specified. Please input a list/tuple of 3 numbers'
#            self.config['rf_atten']=0
#            for att in range(3):
#                set_gain-=self._rf_atten_write(att,0-gain[att])   
#                self.config['rf_attens'][att]=set_gain
#                self.config['rf_atten']+=set_gain
#                if self.config['rf_atten_gain_calfiles'][att] != 'none':
#                    self.config['rf_gain']+=self.get_interpolated_attens(fileName=self.config['rf_atten_gain_calfiles'][att],setpoint=set_gain)
#                else:
#                    self.config['rf_gain']+=set_gain
#
#        elif (gain=='auto') or (gain == 'Auto') or (gain == 'AUTO'):
#            self.auto_gain(print_progress=print_progress)
#
#        else:
#            raise RuntimeError("Could not interpret your gain request. Input a float, a list of 3 values for each attenuator or 'auto'.")
#
#    def _rf_atten_write(self,attenuator,attenuation):
#        """Writes to an RF attenuator. Valid attenuators are 0-2 with a range of 0 to 31.5db of attenuation. Attenuator 0 is closest to the antenna, attenuator 2 is closest to the ADC."""
#        #    'atten0_le_pin' : 5,
#        #    'atten1_le_pin' : 6,
#        #    'atten2_le_pin' : 7,
#        #    'atten_clk_pin' : 3,
#        #    'atten_data_pin': 4}
#        if attenuator == 0:
#            set_att=self.rf_frontend.set_atten_db(le_pin=5,data_pin=4,clk_pin=3,atten_db=attenuation)
#        elif attenuator == 1:
#            set_att=self.rf_frontend.set_atten_db(le_pin=6,data_pin=4,clk_pin=3,atten_db=attenuation)
#        elif attenuator == 2:
#            set_att=self.rf_frontend.set_atten_db(le_pin=7,data_pin=4,clk_pin=3,atten_db=attenuation)
#        else:
#            raise RuntimeError("Invalid attenuator %i. Valid attenuators are %s."%(attenuation,range(0,4)))
#        return set_att

    def adc_amplitudes_get(self):
        """Gets the ADC RMS amplitudes."""
        #TODO: CHECK THESE RETURNS!
        rv = {}
        rv['adc_raw']=self.fpga.read_uint('adc_sum_sq0')
        rv['adc_rms_raw']=numpy.sqrt(rv['adc_raw']/float(self.config['adc_levels_acc_len']))
        rv['adc_rms_mv']=rv['adc_rms_raw']*self.config['adc_v_scale_factor']*1000
        rv['adc_dbm']=ratty2.cal.v_to_dbm(rv['adc_rms_mv']/1000.)
        #backout fe gain
        rv['input_dbm']=rv['adc_dbm']-self.config['rf_gain']
        rv['input_rms_mv']=ratty2.cal.dbm_to_v(rv['input_dbm']*1000)
        return rv

    def status_get(self):
        """Reads and decodes the status register. Resets any error flags after reading."""
        rv={}
        value = self.fpga.read_uint('status0')
        self.ctrl_set(clr_status='pulse')
        return {
                'adc_shutdown':bool(value&(1<<4)),
                'adc_overrange':bool(value&(1<<2)),
                'fft_overrange':bool(value&(1<<1))
                }

    def acc_time_set(self,acc_time=None):
        """Set the accumulation length in seconds. If not specified, use the default from the config file."""
        if acc_time >0:
            self.config['acc_period'] = acc_time
        self.config['n_accs'] = int(self.config['acc_period'] * float(self.config['bandwidth'])/self.config['n_chans'])
        self.logger.info("Setting accumulation time to %2.2f seconds (%i accumulations)."%(self.config['acc_period'],self.config['n_accs']))
        self.fpga.write_int('acc_len',self.config['n_accs'])
        self.ctrl_set(mrst='pulse')

    def acc_time_get(self):
        """Set the accumulation length in seconds"""
        self.config['n_accs'] = self.fpga.read_uint('acc_len')
        self.acc_time=self.config['n_accs']*self.config['n_chans']/float(self.config['bandwidth'])
        self.logger.info("Accumulation time is %2.2f seconds (%i accumulations)."%(self.acc_time,self.config['n_accs']))
        return self.acc_time,self.config['n_accs']

    def get_adc_snapshot(self,trig_level=-1):
        if trig_level>0: 
            self.fpga.write_int('trig_level',trig_level)
            circ_capture=True
        else:
            self.fpga.write_int('trig_level',0)
            circ_capture=False
        return numpy.fromstring(self.fpga.snapshot_get('snap_adc',man_valid=True,man_trig=True,circular_capture=circ_capture,wait_period=-1)['data'],dtype=numpy.int16).byteswap()


    def adc_temp_get(self):
        """Return the temperature in degC for the ADC die."""
        #Note: KATADC and MKADC use same IIC temp sensors; just leverage that here.
        return corr.katadc.get_adc_temp(self.fpga,0)

    def ambient_temp_get(self):
        """Return the temperature in degC for the ambient temperature inside the RATTY2 digital enclosure near the ADC heatsink."""
        #Note: KATADC and MKADC use same IIC temp sensors; just leverage that here.
        return corr.katadc.get_ambient_temp(self.fpga,0)


def ByteToHex( byteStr ):
    """
    Convert a byte string to it's hex string representation e.g. for output.
    """
    
    # Uses list comprehension which is a fractionally faster implementation than
    # the alternative, more readable, implementation below
    #   
    #    hex = []
    #    for aChar in byteStr:
    #        hex.append( "%02X " % ord( aChar ) )
    #
    #    return ''.join( hex ).strip()        

    return ''.join( [ "%02X " % ord( x ) for x in byteStr ] ).strip()


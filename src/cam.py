#!/usr/bin/env python
'''
You need to have KATCP and CORR installed.
Get them from:
http://pypi.python.org/pypi/katcp
and
http://casper.berkeley.edu/svn/trunk/projects/packetized_correlator/corr-0.4.0/

Hard-coded for 32bit unsigned numbers.
\nAuthor: Jason Manley, Feb 2011.
'''

import corr
import time
import numpy
# import struct
import sys
import logging
import ratty2
import cal
# import conf
# import iniparse
# import os
import valon_synth
import socket
# import wx


class spec:
    def __init__(self, connect=False, log_handler=None,
                 log_level=logging.INFO, **kwargs):

        if log_handler is None:
            log_handler = corr.log_handlers.DebugLogHandler(100)

        self.lh = log_handler
        self.logger = logging.getLogger('RATTY2')
        self.logger.setLevel(log_level)
        self.logger.addHandler(self.lh)

        if not 'config_file' in kwargs:
            kwargs['config_file'] = '/etc/ratty2/default_band0'

        self.cal = ratty2.cal.cal(logger=self.logger, **kwargs)
        self.config = self.cal.config

        #skip the first accumulation (0), which contains junk.
        self.last_acc_cnt = 1

        if connect:
            self.connect()

    def connect(self):
        self.logger.info(
            'Trying to connect to ROACH %s on port %i...' %
            (self.config['roach_ip_str'], self.config['katcp_port']))

        self.fpga = corr.katcp_wrapper.FpgaClient(
            self.config['roach_ip_str'],
            self.config['katcp_port'],
            timeout=3,
            logger=self.logger)
        time.sleep(0.2)
        try:
            self.fpga.ping()
            self.logger.info('KATCP connection to FPGA ok.')
        except:
            self.logger.error(
                'KATCP connection failure. Connection to ROACH failed.')
            raise RuntimeError("Connection to FPGA board failed.")

    def auto_gain(self, print_progress=False):
        """
        Try to automatically set the RF attenuators.
        NOT YET IMPLEMENTED FOR RATTY2
        """
        #TODO!!!
        raise RuntimeError('Not yet implemented for RATTY2')

        self.logger.info('Attempting automatic RF gain adjustment...')
        if print_progress:
            print ('Attempting automatic RF gain adjustment...')
        max_n_tries = 10
        n_tries = 0
        tolerance = 1
        rf_gain = self.config['rf_gain_range'][0]
        self.rf_gain_set(rf_gain)

        time.sleep(0.1)

        self.ctrl_set(mrst='pulse',
                      cnt_rst='pulse',
                      clr_status='pulse',
                      flasher_en=True)
        rf_level = self.adc_amplitudes_get()['adc_dbm']

        if self.status_get()['adc_shutdown'] or\
                self.status_get()['adc_overrange']:
            self.logger.error('Your input levels are too high!')
            raise RuntimeError('Your input levels are too high!')

        while (rf_level < self.config[
                'desired_rf_level'
                ] - tolerance or rf_level > self.config[
                'desired_rf_level'
                ] + tolerance) and n_tries < max_n_tries:

            rf_level = self.adc_amplitudes_get()['adc_dbm']
            difference = self.config['desired_rf_level'] - rf_level
            rf_gain = self.rf_status_get()[1] + difference

            log_str =\
                'Gain was %3.1fdB, resulting in an ADC input level of %5.2fdB.\
                    Trying gain of %4.2fdB...' %\
                (self.rf_status_get()[1], rf_level, rf_gain)

            self.logger.info(log_str)
            if print_progress:
                print log_str
            if self.rf_gain < self.config['rf_gain_range'][0]:
                log_str = 'Gain at minimum, %4.2fdB.' %\
                    self.config['rf_gain_range'][0]
                self.logger.warn(log_str)

                if print_progress:
                    print log_str
                self.rf_gain_set(self.config['rf_gain_range'][0])
                break
            elif rf_gain > self.config['rf_gain_range'][1]:
                log_str = 'Gain at maximum, %4.2fdB.' %\
                    self.config['rf_gain_range'][1]
                self.logger.warn(log_str)
                if print_progress:
                    print log_str
                self.rf_gain_set(self.config['rf_gain_range'][1])
                break

            self.rf_gain_set(rf_gain)
            time.sleep(0.1)
            n_tries += 1

        if n_tries >= max_n_tries:
            log_str = 'Auto RF gain adjust failed.'
            self.logger.error(log_str)
            if print_progress:
                print log_str
        else:
            log_str = 'Auto RF gain adjust success.'
            if print_progress:
                print log_str
            self.logger.info(log_str)

    def auto_fft_shift(self):
        """
        Try to automatically set the FFT shift schedule
        """
        self.ctrl_set(mrst='pulse',
                      cnt_rst='pulse',
                      clr_status='pulse',
                      flasher_en=False)
        stat = self.status_get()
        # orig_fft_shift = self.config['fft_shift']
        fft_shift_adj = self.config['fft_shift']
        while not(stat['fft_overrange']):
            fft_shift_adj = fft_shift_adj << 1
            self.fft_shift_set(fft_shift_adj)
            self.ctrl_set(mrst='pulse',
                          cnt_rst='pulse',
                          clr_status='pulse',
                          flasher_en=False)
            stat = self.status_get()

        fft_shift_adj = fft_shift_adj >> 1
        self.fft_shift_set(fft_shift_adj)
        self.config['fft_shift'] = fft_shift_adj & self.config['fft_shift']

        return (fft_shift_adj & self.config['fft_shift'])

    def _rf_band_switch_calc(self, rf_band=None, cal_config=False):
        """
        Calculates the bitmap for the RF switches to select an RF band.
        Select a band between 1 and 4.\n
        1: 0-800 MHz   2: 700 - 1100   2: 1050-1650 MHz  3: 1950 - 2550
        """
        rf_bands = self.config['rf_switch_layout']
        if rf_band is None:
            rf_band = self.config['band_sel']
        assert rf_band <= len(rf_bands), "Requested RF band is out of range" +\
            "(1-%i)" % len(rf_bands)
        bitmap = int(self.config['rf_switch_layout'][rf_band])
        if cal_config:
            bitmap += cal_config
        return rf_band, bitmap

    def fe_set(self, rf_band=None, gain=None, cal_config=False):
        """
        Configures the analogue box:
        selects bandpass filters and adjusts RF attenuators;
        updates global config with changes.\n
        Select a band between 1 and 6 (rough freqs. shown): \n
        \t 1: 0-750 MHz \n
        \t 2: 700-1100 MHz \n
        \t 3: 950-1670 MHz \n
        \t 4: 1560-2050 MHz \n
        \t 5: 1950-2450 MHz \n
        \t 6: 2450-2900 MHz \n
        Valid gain range is -94.5 to 0dB;
        If single value, import atten settings from atten_setting_map \n
        Alternatively, pass a tuple or list to specify the three values
        explicitly. \n
        If no gain is specified, default to whatever's in the config file \n
        Attenuators / switches set in sequence to limit self-generated RFI\n
        """
        self.config['fe_write'] = True  # Flag for front-end ctrl line check
        initial = not(gain)  # TODO Check this is used wisely - maybe just max
        if gain is None:  # set to config_file or existing value
            gain = self.config['rf_atten']
        gain = round(gain * 2.) / 2.
        rf_band_new = rf_band  # record new variables and keep existing config
        attens_interm = self.cal._rf_atten_calc(self.config['rf_atten'])
        att_l = len(attens_interm)
        rf_band, bitmap_interm =\
            self._rf_band_switch_calc(rf_band=self.config['band_sel'])
        attens = self.cal._rf_atten_calc(self.config['max_atten'])  # max attenuation
        bitmap_check = bitmap_interm
        switch_bitmap_offset =\
            numpy.max([int(len(numpy.binary_repr(n)))
                       for n in self.config['rf_switch_layout']])
        for (att, atten) in enumerate(reversed(attens_interm)):
            bitmap_check += (int(-atten * 2)) <<\
                (switch_bitmap_offset + (6 * int(att)))
        for x in range(att_l):  # Set attenuation to maximum, back-to-front
            bitmap = bitmap_interm
            attens_interm[-1 - x] = attens[-1 - x]
            for (att, atten) in enumerate(reversed(attens_interm)):
                bitmap += (int(-atten * 2)) <<\
                    (switch_bitmap_offset + (6 * int(len(attens) - 1 - att)))
            if initial:
                bitmap = 2 ** 32 - 1  # write all lines high first round
                self.fe_write(bitmap)
                bitmap_check = bitmap
                break
            elif bitmap == bitmap_check:
                continue
            else:
                self.fe_write(bitmap)
                bitmap_check = bitmap
        self.cal.update_atten_bandpass(gain=gain)  # updates config['rf_atten']
        rf_band, bitmap_interm =\
            self._rf_band_switch_calc(rf_band=rf_band_new, cal_config=cal_config)
        self.logger.info("Selected RF band %i." % rf_band)
        self.config['band_sel'] = rf_band
        attens_new = self.cal._rf_atten_calc(self.config['rf_atten'])
        for x in range(att_l + 1):  # Toggle switch, then attens (front->back)
            bitmap = bitmap_interm
            if x > 0:
                attens[x - 1] = attens_new[x - 1]
            for (att, atten) in enumerate(attens):
                bitmap += ((int(-atten * 2))) <<\
                    (switch_bitmap_offset + (6 * int(len(attens) - 1 - att)))
            if bitmap == bitmap_check:
                continue
            else:           
                self.fe_write(bitmap)
                bitmap_check = bitmap
        #print '\nBitmap:\t%i\tBitmap Binary:%s\n'\
        #    %(bitmap, numpy.binary_repr(bitmap))


    def fe_write(self, bitmap, read=True):
        '''
        Write bitstream to front-end control board and read back
        register entries, raising RuntimeError if mismatched.
        '''
        self.fpga.write_int('rf_ctrl0', bitmap)
        if read:
            time.sleep(0.4)  # shortest safe delay for succesful readback
            fb_bitmap = self.fpga.read_uint('rf_fb0')
            if bitmap != fb_bitmap:
                err_mssge = 'Front-End Control Error:\n%s\twritten bitstream\
                    \n%s\treturned bitstream' % (bin(bitmap), bin(fb_bitmap))
                self.logger.error(err_mssge)
                print err_mssge
                self.config['fe_write'] = False

    def set_valon(self, freq=None):
        '''
        Configures Valon Synth connected to ROACH2 ser2net gateway on port 7148
        '''

        # The same port as used by the server
        PORT = 7148
        socket._socketobject.read = socket._socketobject.recv
        socket._socketobject.write = socket._socketobject.send
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.settimeout(2)
        s.connect((self.config['roach_ip_str'], PORT))

        valon = valon_synth.Synthesizer(s)
        valon.set_ref_select(0)

        if freq is None:
            freq = self.config['sample_clk']

        #TODO: Check why both synths are set?
        valon.set_frequency(valon_synth.SYNTH_B, freq / 1.e6)
        time.sleep(0.1)
        valon.set_frequency(valon_synth.SYNTH_A, freq / 1.e6)
        time.sleep(0.1)

        freq_chck = valon.get_frequency(valon_synth.SYNTH_B)

        if freq_chck != (freq / 1.e6):
            raise\
                RuntimeError(
                    '\nValon Synthesizer ERROR!\
                    Tried: %4.1f MHz, read back %4.1f MHz' %
                    (freq / 1.e6, freq_chck))
        s.close()

    def initialise(self, skip_program=False, clk_check=False, input_sel='Q',
                   print_progress=False, progressbar=None, outputText=None):
        """Initialises the system to defaults."""
        if not skip_program:
            if print_progress:
                print '\
                    \tConfiguring Valon frequency to %5.1f MHz...'\
                    % (self.config['sample_clk'] / 1e6),

                sys.stdout.flush()

            self.set_valon()

            if print_progress:
                print 'ok'

            # if progressbar is not None:
            #     wx.CallAfter(progressbar.SetGauge, 70)
            #     wx.CallAfter(progressbar.SetText,
            #                  "Configuring Valon and Programming FPGA...")
            #     wx.CallAfter(progressbar.SetGauge, 75)

            if print_progress:
                print '\tProgramming FPGA with %s...' %\
                    (self.config['bitstream']),
                sys.stdout.flush()

            self.fpga.upload_program_bof(
                '/etc/ratty2/boffiles/' + self.config['bitstream'], 3333)

            # if progressbar is not None:
            #     wx.CallAfter(progressbar.SetGauge, 90)

            time.sleep(2)
            if print_progress:
                print 'ok'
            # if progressbar is not None:
            #     wx.CallAfter(progressbar.SetGauge, 100)

        elif print_progress:
            print 'Reprogramming skipped.'

        if clk_check:
            if print_progress:
                print '\tChecking clocks...',
                sys.stdout.flush()
            est_rate = self.clk_check()
            if print_progress:
                print 'ok, %i MHz' % est_rate

        if print_progress:
            print '\tSelecting RF band %i (%i-%i MHz) and...\n\t\
                adjusting attenuators for %4.1fdB total attenuation...' %\
                (self.config['band_sel'],
                 self.config['ignore_low_freq'] / 1.e6,
                 self.config['ignore_high_freq'] / 1.e6,
                 self.config['rf_atten']),
        self.fe_set()
        if print_progress:
            print 'ok'

        if print_progress:
            print '\tTotal frontend gain: %3.1fdb' %\
                ((numpy.mean(self.config['system_bandpass'])))

        if print_progress:
            print '\tConfiguring FFT shift schedule to %i...' %\
                self.config['fft_shift'],
            sys.stdout.flush()

        self.fft_shift_set()
        if print_progress:
            print 'ok'

        if print_progress:
            print '\tConfiguring accumulation period to %4.2f seconds...' %\
                self.config['acc_period'],
            sys.stdout.flush()
        self.acc_time_set(self.config['acc_period'])

        if print_progress:
            print 'ok'

        if print_progress:
            print '\tClearing status and aligning PFB...',
            sys.stdout.flush()

        self.ctrl_set(mrst='pulse',
                      cnt_rst='pulse',
                      clr_status='pulse',
                      flasher_en=False)

        if print_progress:
            print 'ok'

        stat = self.status_get()
        if stat['adc_shutdown']:
            log_msg = 'ADC selfprotect due to overrange!'
            self.logger.error(log_msg)
            if print_progress:
                print log_msg
        elif stat['adc_overrange']:
            log_msg = 'ADC is clipping!'
            self.logger.warn(log_msg)
            if print_progress:
                print log_msg
        elif stat['fft_overrange']:
            log_msg = 'FFT is overflowing!'
            self.logger.error(log_msg)
            if print_progress:
                print log_msg

    def clk_check(self):
        """
        Performs a clock check and returns an estimate of the
        FPGA's clock frequency.
        """
        est_rate = round(self.fpga.est_brd_clk())
        if est_rate > (self.config['fpga_clk'] / 1e6 + 1) or\
                est_rate < (self.config['fpga_clk'] / 1e6 - 1):
            self.logger.error(
                'FPGA clock rate is %i MHz where we expect it to be %i MHz.'
                % (est_rate, self.config['fpga_clk'] / 1e6))

        return est_rate

    def cal_gains(self, low, high):
        """
        Used as part of system calibration, to derive RF
        attenuator calibration files.
        NOT YET IMPLEMENTED FOR RATTY2!
        """

        #TODO Implement for RATTY2
        raise RuntimeError('Not yet implmeneted for RATTY2')
        base_gain = self.rf_atten_set((low + high) / 2.)
        time.sleep(0.2)
        base_power = self.adc_amplitudes_get()['adc_dbm'] - base_gain
        gain_cal = []
        for g in numpy.arange(low,
                              high + self.config['rf_gain_range'][2],
                              self.config['rf_gain_range'][2]):
            self.rf_atten_set(g)
            time.sleep(0.2)
            gain_cal.append(self.adc_amplitudes_get()['adc_dbm'] - base_power)
        return gain_cal

    def get_spectrum(self, last_acc_cnt=None):
        """
        Gets data from ROACH board and returns the spectra and the state
        of the roach at the last timestamp. Units of 'cal_spectrum' are
        dBm unless an antenna was specified in your config file, in which
        case units are dBuV/m.\n
        Performs bandpass correction, fft_scaling adjustment, backs out
        number of accumulations, RF frontend gain etc.\n
        """
        if last_acc_cnt is None:
            last_acc_cnt = self.last_acc_cnt

        while self.fpga.read_uint('acc_cnt') <= (last_acc_cnt):
            #Wait until the next accumulation has been performed.
            time.sleep(0.1)
        spectrum = numpy.zeros(self.config['n_chans'])
        spectrum = self.read_roach_spectrum()

        #for i in range(self.config['n_par_streams']):
        #    spectrum[i::self.config['n_par_streams']] =\
        #        numpy.fromstring(
        #            self.fpga.read('%s%i' % (
        #                self.config['spectrum_bram_out_prefix'], i),
        #                self.config['n_chans'] / self.config['n_par_streams'] * 8),
        #            dtype=numpy.uint64).byteswap()
        #
        #self.last_acc_cnt = self.fpga.read_uint('acc_cnt')
        #if self.config['flip_spectrum']:
        #    spectrum = spectrum[::-1]
        #stat = self.status_get()
        #ampls = self.adc_amplitudes_get()
        #if stat['adc_shutdown']:
        #    self.logger.error('ADC selfprotect due to overrange!')
        #elif stat['adc_overrange']:
        #    self.logger.warning('ADC is clipping!')
        #elif stat['fft_overrange']:
        #    self.logger.error('FFT is overflowing!')

        cal_spectrum = r.cal.calibrate_pfb_spectrum(spectrum)

        return {'raw_spectrum': spectrum,
                'calibrated_spectrum': cal_spectrum,
                'timestamp': time.time(),
                'acc_cnt': self.last_acc_cnt,
                'adc_overrange': stat['adc_overrange'],
                'fft_overrange': stat['fft_overrange'],
                'adc_shutdown': stat['adc_shutdown'],
                'adc_level': ampls['adc_dbm'],
                'input_level': ampls['input_dbm'],
                'adc_temp': self.adc_temp_get(),
                'ambient_temp': self.ambient_temp_get()}

    def fft_shift_set(self, fft_shift_schedule=None):
        """
        Sets the FFT shift schedule (divide-by-two) on each FFT stage.
        Input is an integer representing a binary bitmask for shifting.
        If not specified as a parameter to this function, set the default
        level from the config file.
        """
        if fft_shift_schedule is None:
            fft_shift_schedule = self.config['fft_shift']
        self.fpga.write_int('fft_shift', fft_shift_schedule)
        self.config['fft_shift'] = fft_shift_schedule
        self.config['fft_scale'] = 2 ** (cal.bitcnt(fft_shift_schedule))
        self.logger.info("Set FFT shift to %8x (scaling down by %i)." %
                         (fft_shift_schedule, self.config['fft_scale']))

    def fft_shift_get(self):
        """
        Fetches the current FFT shifting schedule from the hardware.
        """
        self.config['fft_shift'] = self.fpga.read_uint('fft_shift')
        self.config['fft_scale'] = 2 ** (cal.bitcnt(self.config['fft_shift']))
        return self.config['fft_shift']

    def ctrl_get(self):
        """
        Reads and decodes the values from the control register.
        """
        value = self.fpga.read_uint('control')
        return {'mrst': bool(value & (1 << 0)),
                'cnt_rst': bool(value & (1 << 1)),
                'clr_status': bool(value & (1 << 3)),
                'adc_protect_disable': bool(value & (1 << 13)),
                'flasher_en': bool(value & (1 << 12)),
                'raw': value,
                }

    def ctrl_set(self, **kwargs):
        """
        Sets bits of all the Fengine control registers.
        Keeps any previous state.
             \nPossible boolean kwargs:
             \n\t adc_protect_disable
             \n\t flasher_en
             \n\t clr_status
             \n\t mrst
             \n\t cnt_rst
        """

        key_bit_lookup =\
            {'adc_protect_disable': 13,
             'flasher_en': 12,
             'clr_status': 3,
             'cnt_rst': 1,
             'mrst': 0, }

        value = self.ctrl_get()['raw']
        run_cnt = 0
        run_cnt_target = 1
        while run_cnt < run_cnt_target:
            for key in kwargs:
                if (kwargs[key] == 'toggle') and (run_cnt == 0):
                    value = value ^ (1 << (key_bit_lookup[key]))
                elif (kwargs[key] == 'pulse'):
                    run_cnt_target = 3
                    if run_cnt == 0:
                        value = value & ~(1 << (key_bit_lookup[key]))
                    elif run_cnt == 1:
                        value = value | (1 << (key_bit_lookup[key]))
                    elif run_cnt == 2:
                        value = value & ~(1 << (key_bit_lookup[key]))
                elif kwargs[key] is True:
                    value = value | (1 << (key_bit_lookup[key]))
                elif kwargs[key] is False:
                    value = value & ~(1 << (key_bit_lookup[key]))
                else:
                    raise RuntimeError(
                        "Sorry, you must specify True, False,\
                            'toggle' or 'pulse' for %s." % key)
            self.fpga.write_int('control', value)
            run_cnt = run_cnt + 1

    def adc_amplitudes_get(self):
        """
        Gets the ADC RMS amplitudes.
        """
        #TODO: CHECK THESE RETURNS!
        rv = {}
        rv['adc_raw'] = self.fpga.read_uint('adc_sum_sq0')
        rv['adc_rms_raw'] =\
            numpy.sqrt(rv['adc_raw']/float(self.config['adc_levels_acc_len']))
        rv['adc_rms_mv'] =\
            rv['adc_rms_raw']*self.config['adc_v_scale_factor']*1000
        rv['adc_dbm'] = ratty2.cal.v_to_dbm(rv['adc_rms_mv']/1000.)
        rv['input_dbm'] =\
            rv['adc_dbm']-(numpy.mean(self.config['system_bandpass']))

        # CHECK THIS! WHY dBm*1000?
        rv['input_rms_mv'] = ratty2.cal.dbm_to_v(rv['input_dbm']*1000)
        # rv['input_rms_mv']=ratty2.cal.dbm_to_v(rv['input_dbm'])
        return rv

    def status_get(self):
        """
        Reads and decodes the status register.
        Resets any error flags after reading.
        """
        value = self.fpga.read_uint('status0')
        self.ctrl_set(clr_status='pulse')
        return {'adc_shutdown': bool(value & (1 << 4)),
                'adc_overrange': bool(value & (1 << 2)),
                'fft_overrange': bool(value & (1 << 1))}

    def acc_time_set(self, acc_time=None):
        """
        Set the accumulation length in seconds.
        If not specified, use the default from the config file.
        """
        if acc_time > 0:
            self.config['acc_period'] = acc_time

        self.config['n_accs'] =\
            int(self.config['acc_period'] * float(
                self.config['bandwidth'])/self.config['n_chans'])

        self.logger.info(
            "Setting accumulation time to %2.2f seconds (%i accumulations)." %
            (self.config['acc_period'], self.config['n_accs']))

        self.fpga.write_int('acc_len', self.config['n_accs'])
        self.ctrl_set(mrst='pulse')

    def acc_time_get(self):
        """
        Set the accumulation length in seconds
        """
        self.config['n_accs'] = self.fpga.read_uint('acc_len')
        self.acc_time =\
            self.config['n_accs']*self.config['n_chans']/float(
                self.config['bandwidth'])
        self.logger.info(
            "Accumulation time is %2.2f seconds (%i accumulations)." %
            (self.acc_time, self.config['n_accs']))
        return self.acc_time, self.config['n_accs']

    def get_adc_snapshot(self, trig_level=-1, wait_period=-1):
        if trig_level > 0:
            self.fpga.write_int('trig_level', trig_level)
            circ_capture = True
        else:
            self.fpga.write_int('trig_level', 0)
            circ_capture = False

        try:
            result = numpy.fromstring(
                self.fpga.snapshot_get('snap_adc',
                                       man_valid=True,
                                       man_trig=True,
                                       circular_capture=circ_capture,
                                       wait_period=wait_period)['data'],
                                       dtype=numpy.int16).byteswap()
            success = True
        except:
            result = None
            success = False
        return result, success


    def adc_temp_get(self):
        """
        Return the temperature in degC for the ADC
        """
        #Note: KATADC and MKADC use same IIC temp sensors;
        #just leverage that here.
        return corr.katadc.get_adc_temp(self.fpga, 0)

    def ambient_temp_get(self):
        """
        Return the temperature in degC for the ambient temperature
        inside the RATTY2 digital enclosure near the ADC heatsink.
        """
        #Note: KATADC and MKADC use same IIC temp sensors;
        #just leverage that here.
        return corr.katadc.get_ambient_temp(self.fpga, 0)








    def read_roach_spectrum(self):
        spectrum = numpy.zeros(self.config['n_chans'])

        for i in range(self.config['n_par_streams']):
            spectrum[i::self.config['n_par_streams']] =\
                numpy.fromstring(
                    self.fpga.read('%s%i' % (
                        self.config['spectrum_bram_out_prefix'], i),
                        self.config['n_chans'] / self.config['n_par_streams'] * 8),
                    dtype=numpy.uint64).byteswap()
        self.last_acc_cnt = self.fpga.read_uint('acc_cnt')
        if self.config['flip_spectrum']:
            spectrum = spectrum[::-1]
        stat = self.status_get()
        ampls = self.adc_amplitudes_get()
        if stat['adc_shutdown']:
            self.logger.error('ADC selfprotect due to overrange!')
        elif stat['adc_overrange']:
            self.logger.warning('ADC is clipping!')
        elif stat['fft_overrange']:
            self.logger.error('FFT is overflowing!')
        return spectrum



    def noise_calibration_cam(self, last_acc_cnt=None, verbose=True,
                              digital_spectrum=False, plot_cal=False,
                              default_file_path='./'):
        """
        Configures front-end for noise-source input, perform calibration
        measurements (OFF / ON), pass data to cal.py and safe relevant results to file.  
        """
        #print 'rf atten, strt of cam auto cal', self.config['rf_atten']
        noise_cal_nap = self.config['acc_period'] * 2.2
        #print '\nlast acc count:\t', last_acc_cnt
        if last_acc_cnt is None:
            last_acc_cnt = self.last_acc_cnt
        t1 = time.time()
        #print '\nlast acc count:\t', last_acc_cnt
        #print '\ntime 1, last_acc_cnt ', t1, last_acc_cnt
        cal_config = self.config['noise_source_on_layout'] +\
            self.config['input_switch_on_layout']  # ON measurement configuration
        #print '\ncal_config:\t%i\t(bin: %s)\n' %(cal_config, numpy.binary_repr(cal_config))
        self.fe_set(cal_config=cal_config)
        time.sleep(noise_cal_nap)
        # rf_band=self.config['band_sel'], gain=self.config['rf_atten'],
        while self.fpga.read_uint('acc_cnt') <= (last_acc_cnt):
            #Wait until the next accumulation has been performed.
            time.sleep(0.1)
        self.last_acc_cnt = self.fpga.read_uint('acc_cnt')
        print '\n\tself.last_acc_cnt 1) in cam noise cal cam:', self.last_acc_cnt
        hot_spectrum = self.read_roach_spectrum()
        #print '\nhot', time.time(), numpy.mean(hot_spectrum)
        hot_spectrum = self.cal.calibrate_pfb_spectrum(hot_spectrum, noise_cal=True)
        cal_config -= self.config['noise_source_on_layout']  # OFF measurement configuration
        #print '\ncal_config:\t%i\t(bin: %s)\n' %(cal_config, numpy.binary_repr(cal_config))
        #print '\nduration 1, last_acc_cnt ', (time.time()-t1), last_acc_cnt
        t1 = time.time()
        # print "cal config:", cal_config
        self.fe_set(cal_config=cal_config)
        time.sleep(noise_cal_nap)
        #  rf_band=self.config['band_sel'], gain=self.config['rf_atten'],
        #time.sleep(8)
        while self.fpga.read_uint('acc_cnt') <= (last_acc_cnt):
            time.sleep(0.1)  # Wait till next accumulation TODO check continuous impact
        self.last_acc_cnt = self.fpga.read_uint('acc_cnt')
        print '\n\tself.last_acc_cnt 2) in cam noise cal cam:', self.last_acc_cnt
        cold_spectrum = self.read_roach_spectrum()
        cold_spectrum = self.cal.calibrate_pfb_spectrum(cold_spectrum,
                                                        noise_cal=True)
        #print '\nduration 2, last_acc_cnt ', (time.time()-t1), last_acc_cnt
        #print '\ncold (time, mean (dB))', time.time(), numpy.mean(cold_spectrum) 
        #print 'rf atten, before dig spec', self.config['rf_atten']
        if digital_spectrum:
            atten_setting=self.config['rf_atten']
            self.fe_set(gain=self.config['max_atten'],
                        rf_band=self.config['band_sel'])
            #print 'rf atten,  dig spec', self.config['rf_atten']
            time.sleep(noise_cal_nap)
            while self.fpga.read_uint('acc_cnt') <= (last_acc_cnt):
                time.sleep(0.1)  # Wait till next accumulation
            self.last_acc_cnt = self.fpga.read_uint('acc_cnt')
            print '\n\tself.last_acc_cnt 3) in cam noise cal cam:', self.last_acc_cnt
            digital_spectrum = self.read_roach_spectrum()
            digital_spectrum =\
                self.cal.calibrate_pfb_spectrum(digital_spectrum, noise_cal=True)
            #print '\nhot', time.time(), numpy.mean(digital_spectrum)
        #print 'atten setting cam nc', atten_setting
        self.fe_set(gain=atten_setting)
        tsys, gain, tsys_mean, gain_mean =\
            self.cal.noise_calibration_calculation(hot_spectrum, cold_spectrum,
                                                   digital_spectrum=digital_spectrum,
                                                   plot_cal=plot_cal,
                                                   default_file_path=default_file_path)
        #print '\nduration 2, last_acc_cnt ', (time.time()-t1), last_acc_cnt
        #print 'rf atten, end of cam auto cal', self.config['rf_atten']
        if verbose:   
            pnt1 = 6666
            print '\n\nTsys @ %.2f Hz:\t%2f K' %(self.config['freqs'][pnt1],
                                             tsys[pnt1])
            print 'Gain @ %.2f Hz:\t%2f dB' %(self.config['freqs'][pnt1],
                                              10*numpy.log10(gain[pnt1]))
            print "Mean Tsys:\t%.2f K" %tsys_mean
            print "Mean in-band gain:\t%.2f dB\n\n" %(10*numpy.log10(gain_mean))

        self.fe_set(rf_band=self.config['band_sel'], gain=self.config['rf_atten'],
                    cal_config=False)  # return to external rf input TODO CHECK!
        return tsys, gain, tsys_mean, gain_mean 
            #  hot_spectrum, cold_spectrum, digital_spectrum


def ByteToHex(byteStr):
    """
    Convert a byte string to it's hex string representation e.g. for output.
    """
    # Uses list comprehension which is a fractionally faster implementation
    # than the alternative, more readable, implementation below
    #
    #    hex = []
    #    for aChar in byteStr:
    #        hex.append( "%02X " % ord( aChar ) )
    #
    #    return ''.join( hex ).strip()

    return ''.join(["%02X " % ord(x) for x in byteStr]).strip()

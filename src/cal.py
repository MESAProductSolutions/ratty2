# -*- coding: utf-8 -*-
import numpy
import scipy
import scipy.interpolate
# import iniparse
import ratty2
# import corr
import csv

'''
2013-05-03 mods to support programmable attenuator spectral calibration.

2013-04-24 rf_atten and fe_amp rolled into one; arbitrary kwargs to _init now
override config_file params.

smoothing functions from:
http://www.swharden.com/blog/2008-11-17-linear-data-smoothing-in-python/
'''

c = 299792458.  # speed of light in m/s
k = 1.38064852e-23  # Boltzmann constant
# For when you have everything working and ready to install with distutils
# cal_file_path = "/etc/rfi_sys/cal_files/"

# For development
cal_file_path = "/etc/ratty2/cal_files/"


def cal_files(filename):
    if filename[0] == '/' or filename[0] == '.':
        return filename
    else:
        return "%s%s" % (cal_file_path, filename)


def smoothList(list, strippedXs=False, degree=10):
    if strippedXs is True:
        return Xs[0:-(len(list) - (len(list) - degree + 1))]
    smoothed = [0] * (len(list) - degree + 1)
    for i in range(len(smoothed)):
        smoothed[i] = sum(list[i:i+degree])/float(degree)
    return smoothed


def smoothListTriangle(list, strippedXs=False, degree=5):
    weight = []
    window = degree*2-1
    smoothed = [0.0]*(len(list)-window)
    for x in range(1, 2*degree):
        weight.append(degree-abs(degree-x))
    w = numpy.array(weight)
    for i in range(len(smoothed)):
        smoothed[i] = sum(numpy.array(list[i:i+window])*w)/float(sum(w))
    return smoothed


def smoothListGaussian(list, strippedXs=False, degree=5):
    window = degree*2-1
    weight = numpy.array([1.0]*window)
    weightGauss = []
    for i in range(window):
        i = i-degree+1
        frac = i/float(window)
        gauss = 1/(numpy.exp((4*(frac))**2))
        weightGauss.append(gauss)
    weight = numpy.array(weightGauss)*weight
    smoothed = [0.0]*(len(list)-window)
    for i in range(len(smoothed)):
        smoothed[i] = sum(numpy.array(list[i:i+window])*weight)/sum(weight)
    return smoothed


def calc_mkadc_bandpass(chans, adc_cal_file_i):
    gain_mkAdc = getDictFromCSV(adc_cal_file_i)
    x = numpy.arange(0, len(gain_mkAdc), 1)
    y = numpy.zeros(len(gain_mkAdc))
    for i in range(len(gain_mkAdc)-1):
        y[i] = gain_mkAdc[(i)*10.0]
    x_int = numpy.arange(0, len(gain_mkAdc), len(gain_mkAdc)/(chans*1.0))
    bandFlat = numpy.subtract(10, numpy.interp(x_int, x, y))
    return bandFlat


def dmw_per_sq_m_to_dbuv(dbmw):
    # from:
    # http://www.ahsystems.com/notes/RFconversions.php:
    # dBmW/m2 = dBmV/m - 115.8
    return dbmw + 115.8


def dbuv_to_dmw_per_sq_m(dbuv):
    # from http://www.ahsystems.com/notes/RFconversions.php:
    # dBmW/m2 = dBmV/m - 115.8
    return dbuv - 115.8


def dbm_to_dbuv(dbm):
    return dbm + 106.98


def dbuv_to_dbm(dbuv):
    return dbuv - 106.98


def v_to_dbuv(v):
    return 20 * numpy.log10(v * 1e6)


def dbuv_to_v(dbuv):
    return (10 ** (dbuv / 20.)) / 1e6


def dbm_to_v(dbm):
    return numpy.sqrt(10. ** (dbm / 10.) / 1000. * 50.)


def v_to_dbm(v):
    return 10. * numpy.log10(v * v / 50. * 1000.)


def bitcnt(val):
    '''Counts the number of set bits in the binary value.'''
    ret_val = 0
    shift_val = val
    while shift_val >= 1:
        if shift_val & 1:
            ret_val += 1
        shift_val = shift_val >> 1
    return ret_val


def polyfit(freqs, gains, degree=9):
    """ Just calls numpy.polyfit.
        Mostly here as a reminder."""
    return numpy.polyfit(freqs, gains, deg=degree)


def af_from_gain(freqs, gains):
    """Calculate the antenna factor (in dB/m) from a list of frequencies
       (in Hz) and gains (in dBi).

        There are a number of assumptions made for this to be valid:
            1) Far-field only (plane wave excitation).
            2) 50ohm system.
            3) antenna is polarisation matched.
            4) effects of impedance mismatch are included.
    """
    #From Howard's email:
    #The antenna factors are derived from gain by taking:
    #  9.73/(lambda(sqrt(Gain))
    #Note the gain here is the non-dB gain.
    #It is worth noting that in dB’s the conversion is
    #  19.8 – 20log10(lambda) – 10 log(Gain)
    #avoid the divide by zero error:
    if freqs[0] == 0:
        freqs[0] = numpy.finfo(numpy.float).eps
    return 19.8 - 20*numpy.log10(c/freqs) - gains


def gain_from_af(freqs, afs):
    """Calculate the gain (in dBi) from the Antenna Factor (in dB/m)."""
    return 19.8 - 20*numpy.log10(c/freqs) - afs


def getDictFromCSV(filename):
    # import csv
    f = open(filename)
    f.readline()
    reader = csv.reader(f, delimiter=',')
    mydict = dict()
    for row in reader:
        mydict[float(row[0])] = float(row[1])
    return mydict


def get_gains_from_csv(filename):
    freqs = []
    gains = []
    more = True
    fp = open(filename, 'r')
    # import csv
    fc = csv.DictReader(fp)
    while(more):
        try:
            raw_line = fc.next()
            freqs.append(numpy.float(raw_line['freq_hz']))
            gains.append(numpy.float(raw_line['gain_db']))
        except:
            more = False
            break
    return freqs, gains


class cal:
    def __init__(self, logger=None, **kwargs):
        """Either specificy a config_file or else specify individual kwargs."""
        self.config = {}
        if kwargs.has_key('config_file'):
            self.config = ratty2.conf.rattyconf(**kwargs)

        for key in kwargs:
            self.config[key] = kwargs[key]
        self.config['bandwidth'] = self.config['sample_clk']/2
        self.config['chan_width'] =\
            numpy.float(self.config['bandwidth'])/self.config['n_chans']

        self.config['nyquist_zone'] =\
            int(numpy.ceil(
                (self.config['ignore_high_freq'
                             ]+self.config['ignore_low_freq'
                                           ])/2/(self.config['bandwidth'])))

        if bool(self.config['nyquist_zone'] & 1):
            self.config['flip_spectrum'] = False
            self.chan_offset = 0
        else:
            self.config['flip_spectrum'] = True
            self.chan_offset = self.config['bandwidth']/self.config['n_chans']

        # self.start_freq =\  # TODO Remove this chunk?
        #     (self.config['nyquist_zone']-1)*self.config['bandwidth'] +\
        #     chan_offset
        # self.stop_freq =\
        #     (self.config['nyquist_zone'])*self.config['bandwidth'] +\
        #     chan_offset
        # self.config['freqs'] = numpy.linspace(self.start_freq,
        #                                       self.stop_freq,
        #                                       self.config['n_chans'],
        #                                       endpoint=False)

        self.generate_freqs(self.config['n_chans'])

        if (self.config['antenna_bandpass_calfile'] != 'none'):
            self.config['antenna_bandpass'] =\
                self.get_interpolated_gains(
                    self.config['antenna_bandpass_calfile'])
        else:
            self.config['antenna_bandpass'] =\
                numpy.zeros(self.config['n_chans'])

        self.config['ant_factor'] =\
            af_from_gain(self.config['freqs'], self.config['antenna_bandpass'])
        self.config['fft_scale'] = bitcnt(self.config['fft_shift'])

        #self.update_atten_bandpass(gain=self.config['rf_atten'])
        #  Also performed as part of cam init routine 

        #Override any defaults:
        for key in kwargs:
            self.config[key] = kwargs[key]


    def generate_freqs(self, n_chans): 
        start_freq =\
            (self.config['nyquist_zone']-1)*self.config['bandwidth'] +\
            self.chan_offset
        stop_freq =\
            (self.config['nyquist_zone'])*self.config['bandwidth'] +\
            self.chan_offset
        self.config['freqs'] = numpy.linspace(start_freq,
                                              stop_freq,
                                              n_chans,
                                              endpoint=False)

    def update_atten_bandpass(self, gain):
        """
        Extract the RF bandpass calibration from the cal file and
        update global config for the user-specified gain setting.
        Update Ip1dB level to track compression limit
        """
        self.config['rf_atten'] = gain
        try:
            self.config['system_bandpass']
        except:
            self.config['system_bandpass'] = []# initiate list
        if (self.config['system_bandpass_calfile'] != 'none'):
            self.config['system_bandpass'] = self.get_interpolated_attens(
                filename=self.config['system_bandpass_calfile'],
                #PARALLEL ATTENUATORS
                atten_db=self.config['rf_atten'],
                freqs_hz=self.config['freqs'])
        else:
            self.config['system_bandpass'] = numpy.ones(
                self.config['n_chans']) * self.config['rf_atten']
        asm = open(cal_files(self.config['rf_ip1db_map']))
        asmv = csv.DictReader(asm, delimiter=',')
        for line in asmv:
            if float(line['attenuation']) == float(gain):
                rf_ip1db = line['ip1db:limiting_component'].split(':')
                self.config['rf_ip1db'] = [float(rf_ip1db[0]),
                        rf_ip1db[1]]
                break

    def _rf_atten_calc(self, gain=None):
        """Determines the attenuation for each of the 3 RF attenuators
           in the RF box. \n
        \t Valid gain reange is config['max_atten'] to 0dB;
           Single gain value settings read from atten_setting_map\n
        \t Alternatively, pass a tuple or list to specify the three values
           explicitly. \n
        \t If no gain is specified, default to whatever's in the config
           file \n"""
        if type(gain) == list or\
                type(gain) == numpy.ndarray or\
                type(gain) == tuple:
            ''' Make rf_attens equal to specified gain '''
            rf_attens = gain

        elif type(gain) == int or type(gain) == float:
            asm = open(cal_files(self.config['atten_setting_map']))
            asmv = csv.DictReader(asm, delimiter=',')
            gain = round(gain*2)/2
            for line in asmv:
                if float(line['Attenuation']) == float(gain):
                    rf_attens_temp =\
                        line['Att1_Sttng:Att2_Sttng:Att3_Sttng'].split(':')
                    rf_attens = [float(rf_attens_temp[0]),
                                 float(rf_attens_temp[1]),
                                 float(rf_attens_temp[2])]
                    break
            asm.close()
        elif gain is None:
            if self.config.has_key('max_atten'):  # Set to max attenuation
                rf_attens = [self.config['max_atten']/3. for att in range(3)]
            #if self.config.has_key('rf_atten'):  # else set to previous setting
            #    rf_attens = [self.config['rf_atten']/3. for att in range(3)]
            #elif self.config.has_key('rf_attens'):
            #    rf_attens = self.config['rf_attens']
            else:
                raise RuntimeError(
                    'Unable to figure out your config file\'s frontend gains;')
        else:
            raise RuntimeError(
                'Unable to figure out your requested attenuation;')
        assert len(rf_attens) == 3, 'Incorect number of gains specified.\
            Please input a list/tuple of values equal to the no. of attens'
        for att in range(3):
            assert (-31.5 <= rf_attens[att] <= 0),\
                "Range for attenuator %i (%3.1f) out of range (-31.5 to 0)." %\
                (att, rf_attens[att])
        assert sum(rf_attens) >= self.config['max_atten'], '\nRequested' +\
            'attenuation out of range.\nLimit attenuation to max:\t' +\
            str(self.config['max_atten']) + ' dB'
        #PARALLEL ATTENUATORS ONLY STEP IN 1dB STEPS - NEED TO ROUND 1dB
        return [round(att*2)/2 for att in rf_attens]

    def plot_bandshape(self, freqs):
        import pylab
        #TODO: Undefined "bandshape"
        pylab.plot(self.config['system_bandpass'](freqs))
        pylab.title('Bandpass calibration profile')
        pylab.xlabel('Frequency (Hz)')
        pylab.ylabel('Relative response (dB)')

    def plot_atten_gain_map(self):
        import pylab
        #TODO: Undefined "atten_gain_map"
        inputs = atten_gain_map.keys()
        inputs.sort()
        pylab.plot(inputs, [atten_gain_map[k] for k in inputs])
        pylab.title('RF attenuator mapping')
        pylab.xlabel('Requested value (dB)')
        pylab.ylabel('Actual value (dB)')
  
    def get_interpolated_gains(self, filename):
        """
        Retrieves antenna gain mapping from bandpass csv files and
        interpolates data to return values at 'freqs'.
        """
        cal_freqs, cal_gains = get_gains_from_csv(cal_files(filename))
        inter_freqs = scipy.interpolate.interp1d(cal_freqs,
                                                 cal_gains,
                                                 kind='linear')
 
        return inter_freqs(self.config['freqs'])

    def get_interpolated_attens(self, filename, atten_db, freqs_hz):
        freqs = []
        gains = []
        more = True
        fp = open(cal_files(filename), 'r')
        fc = csv.DictReader(fp, delimiter=',')
        while(more):
            try:
                raw_line = fc.next()
            except:
                more = False
                break
            if float(atten_db) == float(raw_line['setpoint_db']):
                for pair in raw_line['freq_hz:measured_gain_db'].split(';'):
                    freq, gain = pair.split(':')
                    freqs.append(float(freq))
                    gains.append(float(gain))
        if not(freqs[0] <= freqs_hz[0]):
            freqs.insert(0, 0.0)
            gains.insert(0, gains[0])
        if not(freqs[len(freqs)-1] >= freqs_hz[len(freqs_hz)-1]):
            freqs.append(freqs_hz[len(freqs_hz)-1])
            gains.append(gains[len(gains)-1])
        inter_attens = scipy.interpolate.interp1d(freqs,
                                                  gains,
                                                  kind='linear')

        return inter_attens(freqs_hz).reshape(len(freqs_hz))

    def plot_ant_gain(self):
        """Plots the antenna gain."""
        import pylab
        pylab.plot(self.config['freqs']/1e6, self.config['antenna_bandpass'])
        pylab.title('Antenna gain %s' %
                    self.config['antenna_bandpass_calfile'])
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Relative response (dBi)')

    def plot_ant_factor(self):
        """
        Plots the antenna factor over the given frequencies as
        calculated from the specified antenna CSV file.
        """
        import pylab
        pylab.plot(self.config['freqs']/1e6, self.config['ant_factor'])
        pylab.title(
            'Antenna factor as a function of frequency (%s)' % self.config[
                'antenna_bandpass_calfile'])
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Antenna factor (dBuV/m)')

    def freq_to_chan(self, frequency, n_chans=None):
        """
        Returns the channel number where a given frequency is to be found.
        Frequency is in Hz.
        """
        #if frequency<0:
        #    frequency=self.config['bandwidth']+frequency
        return\
            int(round(float(frequency) / self.config['bandwidth'] *
                self.config['n_chans']) % self.config['n_chans'])

    def get_input_scale_factor(self):
        """
        Provide the calibration factor to get from an ADC input
        voltage to the actual frontend input voltage. Does not
        perform any frequency-dependent calibration.
        """
        x_low = x_high = None
        for x in range(len(self.config['freqs'])):
            if float(self.config['freqs'][x]) >\
                float(self.config['ignore_low_freq']) and\
                    (x_low is None):
                x_low = x
            elif float(self.config['freqs'][x]) >\
                float(self.config['ignore_high_freq']) or\
                (x == (len(self.config['freqs'])-1)):
                x_high = x
                passband_mean =\
                    numpy.mean(self.config['system_bandpass'][x_low:x_high])
                break

        return 1./(10.**(passband_mean/20.))

    def calibrate_adc_snapshot(self, raw_data):
        """
        Calibrates a raw ADC count timedomain snapshot.
        Returns ADC samples in V, ADC spectrum in dBm,
        input spectrum in dBm and input spectrum of
        n_chans in dBm.
        """
        ret = {}
        ret['adc_mixed'] = numpy.array(raw_data)
        if self.config['flip_spectrum']:
            ret['adc_mixed'][::2] *= -1
        ret['adc_v'] = ret['adc_mixed']*self.config['adc_v_scale_factor']
        ret['input_v'] = ret['adc_v']*self.get_input_scale_factor()
        n_accs = len(raw_data)/self.config['n_chans']/2
        window = numpy.hamming(self.config['n_chans']*2)
        spectrum = numpy.zeros(self.config['n_chans'])
        ret['n_accs'] = n_accs
        # if n_accs < 1:
        #     n_accs = 1
        for acc in range(n_accs):
            # print "for acc ", acc
            spectrum +=\
                numpy.abs((numpy.fft.rfft(
                    ret['adc_v'][self.config['n_chans']*2*acc:\
                                 self.config['n_chans']*2*(acc+1)]*window)[0:self.config['n_chans']]))
        ret['adc_spectrum_dbm'] =\
            20*numpy.log10(spectrum/n_accs/self.config['n_chans']*6.14)
        ret['input_spectrum_dbm'] = ret[
            'adc_spectrum_dbm']-(self.config['system_bandpass'])
        if self.config['antenna_bandpass_calfile'] != 'none':
            ret['input_spectrum_dbuv'] = dbm_to_dbuv(ret[
                'input_spectrum_dbm']) + self.config['antenna_factor']
        return ret







    def omit_bins(self, data, omit):
        try:
            data = [e for i, e in enumerate(data) if i not in omit]
            return numpy.array(data)
        except TypeError:
            return data


    def find_ignore_freq_bins(self):
        '''
        Determine bin numbers for the lower and upper frequency points to ignore 
        '''
        x_low = x_high = None
        for x in range(len(self.config['freqs'])):
            if float(self.config['freqs'][x]) >\
                float(self.config['ignore_low_freq']) and\
                    (x_low is None):
                x_low = x
            elif (float(self.config['freqs'][x]) >\
                float(self.config['ignore_high_freq']) or\
                (x == (len(self.config['freqs'])-1))) and\
                x_high is None:
                x_high = x
                break
        return x_low, x_high

    
    def calibrate_pfb_spectrum(self, spectrum, noise_cal=False):
        '''
        Apply calibration factors to RTA spectrum
        '''
        spectrum /= float(self.config['n_accs'])
        spectrum *= self.config['fft_scale']
        spectrum *= self.config['adc_v_scale_factor']
        spectrum = 10 * numpy.log10(spectrum)
        spectrum -= self.config['pfb_scale_factor']
        if not(noise_cal):
            spectrum -= self.config['system_bandpass']
            if self.config['antenna_bandpass_calfile'] != 'none':
                spectrum = ratty2.cal.dbm_to_dbuv(spectrum)
                spectrum += self.config['ant_factor']
        return spectrum


    def noise_calibration_calculation(self, hot_spectrum, cold_spectrum, 
                                      digital_spectrum = False, plot_cal=False,
                                      default_file_path='./'):
        '''
        Calculation of receiver parameters based on noise-source measurements:
        Y = Non / Noff = hot_spectrum / cold_spectrum
        Tsys = (Tson + Y*Tsoff) / (Y -1)
        G = Pon / (k * B * (Tson + Trx))
        Tson = 10**(NF/10)*T0 + Tamb  # On temp of source
        '''
        t0 = 290.  # Reference temperature
        tamb = 290.  # ambient temperature
        noise_source_NF =\
            self.get_interpolated_gains(
                cal_files(self.config['noise_source_calfile']))
        tson = 10. ** (self.omit_bins(noise_source_NF,
                                      self.config['self_RFI_bins'])
                       / 10.) * t0 + tamb  # Checked 
        y_fact = 10.**(hot_spectrum / 10.) / 10.**(cold_spectrum/10.)
        tsys = (tson - y_fact * tamb) / (y_fact - 1.)  # checked
        ton = 10**((hot_spectrum -
                    self.omit_bins(self.config['system_bandpass'],
                                   r.config['self_RFI_bins'])
                    - 30.) / 10.)\
            / (k * self.config['chan_width'])           
        toff = 10**((cold_spectrum -
                     self.omit_bins(self.config['system_bandpass'],
                                    self.config['self_RFI_bins'])
                     - 30.) / 10.)\
            / (k * self.config['chan_width'])
        y_fact_check = ton / toff  # alternative calculation as check
        gain = 10.**((hot_spectrum - 30.) / 10.)\
            / (k * self.config['chan_width'] * (tsys + tson))  # checked
        low_bin, high_bin = self.find_ignore_freq_bins()
        diff_spectrum = hot_spectrum - cold_spectrum
        mean_diff = numpy.mean(diff_spectrum[low_bin:high_bin])
        tsys_mean = numpy.mean(tsys[low_bin:high_bin])
        gain_mean = numpy.mean(gain[low_bin:high_bin])
        if plot_cal:
            self.plot_noise_calibration(tsys, gain, low_bin, high_bin,
                                        cold_spectrum, hot_spectrum,
                                        digital_spectrum, mean_diff,
                                        default_file_path=default_file_path)
        return tsys, gain, tsys_mean, gain_mean


    def plot_noise_calibration(self, tsys, gain, low_bin, high_bin,
                               cold_spectrum, hot_spectrum,
                               digital_spectrum, mean_diff, save_fig=True,
                               default_file_path='./', diff=True):
        """Plot individual cal. measurement results and spectra."""
        import matplotlib.pyplot as plt
        if not(diff):
            f, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(9, 9))
        else:
            f, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(9, 9))
        ax1.plot(self.config['freqs'][low_bin:high_bin]/1.e6,
                 tsys[low_bin:high_bin])
        ax1.set_xlabel('Frequency (MHz)')
        ax1.set_ylabel('Receiver temperature (K)')
        ax2.plot(self.config['freqs'][low_bin:high_bin]/1.e6,
                 10*numpy.log10(gain[low_bin:high_bin]))
        ax1.grid()
        ax2.set_ylabel('Gain (dB)')
        ax2.set_xlabel('Frequency (MHz)')
        ax2.set_ylabel('Gain (dB)')
        ax2.grid()
        ax3.plot(self.config['freqs'][low_bin:high_bin]/1.e6,
                 hot_spectrum[low_bin:high_bin], label='Hot')
        ax3.plot(self.config['freqs'][low_bin:high_bin]/1.e6,
                 cold_spectrum[low_bin:high_bin], label='Cold')
        dig_title_string = ''  # string to add to title
        if type(digital_spectrum) != bool:
            ax3.plot(self.config['freqs'][low_bin:high_bin]/1.e6,
                     digital_spectrum[low_bin:high_bin],
                     label='Digital Noise Floor')
            mean_dig_cold_diff = numpy.mean(cold_spectrum[low_bin:high_bin] -
                                            digital_spectrum[low_bin:high_bin])
            dig_title_string += '%.2fdB mean Cold-Digital Diff'\
                % mean_dig_cold_diff
        else:
            dig_title_string = ''
        ax3.set_ylabel('Power Spectrum (dBm)')
        ax3.set_xlabel('Frequency (MHz)')
        ax3.legend()
        ax3.grid()
        if diff:
            ax4.plot
            ax4.plot(self.config['freqs'][low_bin:high_bin]/1.e6,
                     hot_spectrum[low_bin:high_bin] - 
                     cold_spectrum[low_bin:high_bin],
                     label='Hot-Cold Difference')
            ax4.set_ylabel('Hot-Cold Difference (dB)')
            ax4.set_xlabel('Frequency (MHz)')
            ax4.legend()
            ax4.grid()

        f.suptitle('Auto-Calibration Spectra and Results\n\n' +\
                   'Band %i, %.1f dB Attenuation Setting\n'\
                   % (self.config['band_sel']+1, self.config['rf_atten']) +\
                   ' (%.2fs Acc. period  ||  %.2fdB mean Hot-Cold Diff.  ||  '\
                   % (self.config['acc_period'], mean_diff) +\
                   dig_title_string + ')')
        if save_fig:
            plt.savefig(default_file_path + 'band' + str(self.config['band_sel'] + 1) +
                        '_min' + str(numpy.abs(self.config['rf_atten'])).split('.')[0] +
                        'dB_atten_auto_cal_spectra_and_results.png',
                        format="png", dpi=600)
        plt.show()
        plt.close()
        return

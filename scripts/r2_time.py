#!/usr/bin/env python2.7

'''
Plots a histogram and time-domain sample of the ADC values
from a specified antenna and pol.

'''
'''
Revisions:
2012-06-15  JRM Update to use objectified package
                trig scale factor now referenced to input levels.
                Added option to not plot histogram.
2011-03-xx  JRM Misc modifications, feature additions etc
2011-02-24  JRM Port to RFI system
2010-12-11: JRM Add printout of number of bits toggling in ADC.
                Add warning for non-8bit ADCs.
2010-08-05: JRM Mods to support variable snap block length.
1.1 PVP Initial.\n

'''

#TODO: Add duty-cycle measurement support.
#TODO: Add trigger count support.

import matplotlib
matplotlib.use('TkAgg')
import ratty2
import time
# import corr
import numpy
# import struct
import sys
import logging
import pylab
import h5py
# import os
# import iniparse
import csv
import ast

global cnt
cnt = 0


def exit_fail():
    print 'FAILURE DETECTED. Log entries:\n',
    try:
        f.flush()
        f.close()
        r.lh.printMessages()
        r.fpga.stop()
    except:
        pass
    if verbose:
        raise
    exit()


def exit_clean():
    try:
        print "Closing file."
        f.flush()
        f.close()
        r.fpga.stop()
    except:
        pass
    exit()


# callback function to draw the data for all the required polarisations
def drawDataCallback(n_samples, trig_level):
    dat = getUnpackedData()
    dat.update(co.calibrate_adc_snapshot(raw_data=dat['adc_raw']))
    calData = dat['input_v']*1000  # in mV
    #freqs=dat['freqs']
    abs_levs = numpy.abs(calData)
    max_lev = numpy.max(abs_levs)
    trigs = numpy.ma.flatnotmasked_edges(
        numpy.ma.masked_less_equal(abs_levs,
                                   (trig_level-4) *
                                   trig_scale_factor))
    if (trigs is None or trigs[0] == 0) and\
       (trig_level > 0) and\
       (max_lev/trig_scale_factor) < trig_level:
        max_pos = numpy.argmax(calData)
        print('ERROR: we asked for a trigger level of %4.2f mV and\
               the hardware reported success, but the maximum level in\
               the returned data was only %4.2fmV.'
              % (trig_level *
                  trig_scale_factor, max_lev))

    if trigs is None:
        max_pos = numpy.argmax(calData)
    else:
        max_pos = trigs[0]

    next_subplot = 0
    if opts.plot_hist:
        subplots[0].cla()
        subplots[0].set_xticks(range(-500, 501, 100))
        histData, bins, patches = subplots[0].hist(dat['adc_raw'],
                                                   bins=256*4,
                                                   range=(-128*4, 128*4))

        if dat['adc_overrange'] or dat['adc_shutdown']:
            subplots[0].set_title('Histogram as at %s'
                                  % (time.ctime(dat['timestamp'])),
                                  bbox=dict(facecolor='red', alpha=0.5))
        else:
            subplots[0].set_title('Histogram as at %s'
                                  % (time.ctime(dat['timestamp'])))
        subplots[0].set_ylabel('Counts')
        subplots[0].set_xlabel('ADC sample bins.')
        matplotlib.pyplot.ylim(ymax=(max(histData) * 1.05))
        next_subplot += 1
    subplots[next_subplot].cla()
    t_start = max(0, max_pos-n_samples/2)
    t_stop = min(len(calData), max_pos+n_samples/2)
    p_data = calData[t_start:t_stop]
    x_range = numpy.arange(t_start-max_pos, t_stop-max_pos)*1.e9/sample_clk

    subplots[next_subplot].plot(x_range, p_data)
    subplots[next_subplot].set_xlim(-n_samples/2*1.e9/sample_clk,
                                    n_samples/2*1.e9/sample_clk)

    if dat['adc_overrange'] or dat['adc_shutdown']:
        subplots[next_subplot].set_title('Time-domain [%i] (max >%4.2fmV)'
                                         % (cnt-1, max_lev),
                                         bbox=dict(facecolor='red', alpha=0.5))
    else:
        subplots[next_subplot].set_title(
            'Time-domain [%i] (max %4.2fmV; ADC %i)'
            % (cnt-1, max_lev, numpy.max(numpy.abs(dat['adc_raw']))))

    subplots[next_subplot].set_ylim(-max_lev-1, max_lev+1)
    subplots[next_subplot].set_ylabel('mV')
    subplots[next_subplot].set_xlabel('Time (nanoseconds).')
    next_subplot += 1

    subplots[next_subplot].cla()
    t_start = 0

    if co.config['antenna_bandpass_calfile'] == 'none':
        calSpectrum = dat['input_spectrum_dbm']
        emptySpectrum = co.calibrate_adc_snapshot(
            raw_data=dat['adc_raw'][t_start:max_pos-1])['input_spectrum_dbm']
    else:
        calSpectrum = dat['input_spectrum_dbuv']
        emptySpectrum = co.calibrate_adc_snapshot(
            raw_data=dat['adc_raw'][t_start:max_pos-1])['input_spectrum_dbuv']

    pylab.hold(True)
    subplots[next_subplot].plot(
        freqs[chan_low:chan_high]/1e6,
        calSpectrum[chan_low:chan_high],
        label='Signal on')
    pylab.hold(True)
    subplots[next_subplot].plot(
        freqs[chan_low:chan_high]/1e6,
        emptySpectrum[chan_low:chan_high],
        label='Quiescent')
    subplots[next_subplot].legend()
    subplots[next_subplot].set_title('Spectrum of capture (%i samples)'
                                     % (len(dat[
                                         'adc_raw'][chan_low:chan_high])))
    subplots[next_subplot].set_ylabel('Level (%s)' % units)
    subplots[next_subplot].set_xlabel('Frequency (MHz)')

    if opts.csv_file:
        csv_writer(dat, quiescent=emptySpectrum)

    fig.canvas.draw()
    if wait_keypress:
        print '\t Press enter to get another capture...'
        raw_input()
    fig.canvas.manager.window.after(100,
                                    drawDataCallback,
                                    n_samples,
                                    trig_level)


def getUnpackedData():
    global cnt
    if play_filename is None:
        print '\nFetching data from ROACH...', co.config['trig_level']
        unpackedBytes = r.get_adc_snapshot(trig_level=co.config['trig_level'])
        print 'done'
        stat = r.status_get()
        ampls = r.adc_amplitudes_get()
        stat['adc_raw'] = unpackedBytes
        stat['adc_level'] = ampls['adc_dbm']
        stat['input_level'] = ampls['input_dbm']
        stat['adc_temp'] = r.adc_temp_get()
        stat['ambient_temp'] = r.ambient_temp_get()
        stat['file_cnt'] = cnt
        stat['timestamp'] = time.time()
        cnt = filewrite(stat)
    else:
        if cnt+1 >= f['adc_raw'].shape[0]:
            print 'No more data; end of file... bye!'
            exit_clean()
        stat = {'adc_overrange': f['adc_overrange'][cnt],
                'adc_raw': f['adc_raw'][cnt],
                'adc_shutdown': f['adc_shutdown'][cnt],
                'adc_level': f['adc_level'][cnt],
                'input_level': f['input_level'][cnt],
                'timestamp': f['timestamp'][cnt],
                'file_cnt': f['file_cnt'][cnt],
                'adc_temp': f['adc_temp'][cnt]}
        cnt += 1

    print '[%i] %s: input level: %5.2f dBm (ADC %5.2f dBm), %f degC.' %\
        (cnt-1, time.ctime(stat['timestamp']), stat['input_level'],
         stat['adc_level'], stat['adc_temp']),
    if stat['adc_shutdown']:
        print 'ADC selfprotect due to overrange!',
    elif stat['adc_overrange']:
        print 'ADC is clipping!',
    else:
        print 'all ok.',
    print ''
    return stat


def filewrite(stat):
    cnt = f['adc_raw'].shape[0]
    print 'Storing entry %i...' % (cnt-1),
    sys.stdout.flush()
    #return #skip for now!
    for ky in stat:
        try:
            #print 'Trying to store %s: '%stat,status[stat]
            f[ky].resize(cnt+1, axis=0)
            f[ky][cnt-1] = stat[ky]
        except KeyError:
            #print 'Creating dataset to store %s: '%ky,stat[ky]
            f.create_dataset(ky, shape=[1], maxshape=[None], data=stat[ky])
    print 'done'
    return cnt


def csv_writer(dat, quiescent):
    #TODO: FIX THIS FOR RATTY2
    raise RuntimeError('Sorry, CSV export not yet implemented for RATTY2')
    fcp = open(str(timestamp)+'.csv', 'w')
    fc = csv.writer(fcp)
    for key in r.config.config.keys():
        if type(r.config[key]) == list:
            fc.writerow([key] + r.config[key])
        elif type(r.config[key]) == numpy.ndarray:
            fc.writerow([key] + r.config[key].tolist())
        else:
            fc.writerow([key] + [r.config[key]])
    fc.writerow(['trig_level']+[trig_level])
    fc.writerow(['timestamp']+['%s' % time.ctime(timestamp)])
    fc.writerow(['adc_overrange'] + [status['adc_overrange']])
    fc.writerow(['fft_overrange'] + [status['fft_overrange']])
    fc.writerow(['adc_shutdown'] + [status['adc_shutdown']])
    fc.writerow(['ave_adc_level_dbm'] + [status['adc_level']])
    fc.writerow(['ave_input_level_dbm'] + [status['input_level']])

    fc.writerow(['raw_adc', 'adc_v', 'input_v',
                 'freq', 'input_spectrum_dbm',
                 'input_spectrum_dbuv', 'quiescent'])

    for i in range(len(cald['adc_v'])):
        if i < n_chans:
            if co.config['antenna_bandpass_calfile'] != 'none':
                fc.writerow(
                    [cald['adc_raw'][i],
                     cald['adc_v'][i],
                     cald['input_v'][i],
                     cald['freqs'][i],
                     cald['input_spectrum_dbm'][i],
                     cald['input_spectrum_dbuv'][i],
                     quiescent[i]])
            else:
                fc.writerow(
                    [cald['adc_raw'][i],
                     cald['adc_v'][i],
                     cald['input_v'][i],
                     cald['freqs'][i],
                     cald['input_spectrum_dbm'][i],
                     quiescent[i]])
        else:
            fc.writerow(
                [cald['adc_raw'][i],
                 cald['adc_v'][i],
                 cald['input_v'][i]])
    fcp.close()


def parseargs(args):
    ret = {}
    for arg in args:
        arg = arg.split('=')
        try:
            ret[arg[0]] = ast.literal_eval(arg[1])
        except ValueError:
            ret[arg[0]] = arg[1]
        except SyntaxError:
            ret[arg[0]] = arg[1]
    return ret


if __name__ == '__main__':
    from optparse import OptionParser
    p = OptionParser()
    p.set_usage('%prog [options] LOG_MESSAGE')
    p.add_option('-v', '--verbose', dest='verbose',
                 action='store_true', default=False,
                 help='Enable debug mode.')
    p.add_option('-o', '--plot_hist', dest='plot_hist',
                 action='store_false', default=True,
                 help='Do not plot the histogram.')
    p.add_option('-t', '--capture_len', dest='capture_len',
                 type='int', default=100,
                 help='Plot this many nano-seconds around the trigger point.\
                       Default:100')
    p.add_option('-n', '--n_chans', dest='n_chans', type='int',
                 default=512,
                 help='Number of frequency channels to resolve in software FFT\
                       Default:512')
    p.add_option('-l', '--trig_level', dest='trig_level', type='float',
                 default=0.,
                 help='Ask the hardware to wait for a signal with at least\
                       this amplitude in mV before capturing.\
                       Default: 0 (disabled, just plot current input).')
    p.add_option('-f', '--play_file', dest='play_file', type='string',
                 default=None,
                 help='Open an existing file for analysis.')
    p.add_option('-s', '--csv_file', dest='csv_file', action='store_true',
                 default=False,
                 help='Output (convert) each timestamp to a separate CSV file.\
                      ')
    p.add_option('-w', '--wait_keypress', dest='wait_keypress',
                 action='store_true', default=False,
                 help='Wait for a user keypress before storing/plotting the\
                       next update.')
    p.add_option('-p', '--skip_prog', dest='fpga_prog', action='store_false',
                 default=True,
                 help='Skip reprogramming the FPGA.')

    p.set_description(__doc__)
    opts, args = p.parse_args(sys.argv[1:])
    kwargs = parseargs(args)

    verbose = opts.verbose
    play_filename = opts.play_file
    wait_keypress = opts.wait_keypress


try:
    kwargs['n_chans'] = opts.n_chans
    if play_filename is None:
        r = ratty2.cam.spec(**kwargs)
        co = r.cal
        print 'Config file %s parsed ok!' % (r.config['config_file'])
        print 'Connecting to ROACH %s...' % r.config['roach_ip_str'],
        r.connect()

        if verbose:
            r.logger.setLevel(logging.DEBUG)
        else:
            r.logger.setLevel(logging.INFO)
        print 'done.'

        r.initialise(skip_program=(not opts.fpga_prog), print_progress=True)
        co.config['trig_scale_factor'] =\
            co.get_input_scale_factor()*co.config['adc_v_scale_factor']

        co.config['trig_level'] =\
            int((opts.trig_level/1000.)/co.config['trig_scale_factor'])
        co.config['n_samples'] =\
            int(opts.capture_len/1.e9*co.config['sample_clk'])

        filename = str(int(time.time())) + ".time.h5"
        print 'Starting file %s.' % filename
        f = h5py.File(filename, mode='w')

        print 'fetching baseline...',
        sys.stdout.flush()
        baseline = r.get_adc_snapshot()
        print 'done'

        usrlog = ('Starting file at %s (%i).'
                  % (time.ctime(), int(time.time())))

        f['/'].attrs['usrlog'] = usrlog

        co.config['usrlog'] = ('Starting file at %i.'
                               % (int(time.time()))).join(args)
        f.create_dataset('adc_raw',
                         shape=[1,
                                len(baseline)],
                         dtype=numpy.int16,
                         maxshape=[None, len(baseline)])

        for key in r.config.config.keys():
            try:
                f['/'].attrs[key] = r.config[key]
            except:
                try:
                    f[key] = r.config[key]
                except TypeError:
                    if r.config[key] is None:
                        f['/'].attrs[key] = 'none'
                    elif type(r.config[key]) == dict:
                        f[key] = r.config[key].items()
    else:
        print 'Opening file %s...' % play_filename
        f = h5py.File(play_filename, 'r')
        conf_ovr = dict(f['/'].attrs)
        for key in f.keys():
            if not key in ['adc_raw',
                           'timestamp',
                           'adc_overrange',
                           'fft_overrange',
                           'adc_shutdown',
                           'adc_level',
                           'input_level']:

                print 'trying', key
                if len(f[key]) > 1:
                    conf_ovr[key] = f[key][:]
                else:
                    conf_ovr[key] = f[key]
        co = ratty2.cal.cal(**conf_ovr)

    freqs = co.config['freqs']
    n_samples = co.config['n_samples']
    sample_clk = co.config['sample_clk']
    trig_scale_factor = co.config['trig_scale_factor']
    trig_level = co.config['trig_level']
    n_chans = co.config['n_chans']

    if co.config['antenna_bandpass_calfile'] == 'none':
        units = 'dBm'
    else:
        units = 'dBuV/m'

    print 'USRLOG: %s' % co.config['usrlog']
    print 'Triggering at a level of %4.2fmV (ADC level of %i).'\
        % (trig_level*trig_scale_factor, trig_level)
    print 'Plotting %i samples.' % n_samples
    chan_low = co.freq_to_chan(co.config['ignore_low_freq'], n_chans=n_chans)
    chan_high = co.freq_to_chan(co.config['ignore_high_freq'], n_chans=n_chans)
    print 'Working with channels %i (%5.1fMHz) to %i (%5.1fMHz).'\
        % (chan_low, freqs[chan_low]/1.e6,
           chan_high, freqs[chan_high]/1.e6)

    # create the subplots
    fig = matplotlib.pyplot.figure()
    subplots = []

    n_subplots = 2
    if opts.plot_hist:
        n_subplots += 1

    for p in range(n_subplots):
        subPlot = fig.add_subplot(n_subplots, 1, p + 1)
        subplots.append(subPlot)

    # start the process
    print 'Starting plots...'
    fig.subplots_adjust(hspace=0.8)
    fig.canvas.manager.window.after(100, drawDataCallback, n_samples,
                                    trig_level)
    matplotlib.pyplot.show()

except KeyboardInterrupt:
    exit_clean()
except:
#    exit_fail()
    raise
print 'Done with all.'
exit_clean()

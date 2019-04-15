#!/usr/bin/env python

'''
Writes spectrum and GPS data from an RFI monitoring spectrometer.\n

'''
#Revisions:\n
#2019-04-15  AJO GPS Threading
#2017-03-24  AJO GPS Coordinates and Data Dump to RPi2
#2012-06-14  JRM Overhaul to more object-oriented code.
#                Diff plots now have maxhold y-axis and peak annotations.
#                Added option to only plot one capture (useful for interacting with a single spectrum)
#2011-03-17  JRM Removed lin plots, added cal
#2011-03-xx  JRM Added various features - logging to file, lin/log plots, status reporting etc.
#2011-02-24  JRM Port to RFI system
#2010-12-11: JRM Add printout of number of bits toggling in ADC.
#                Add warning for non-8bit ADCs.
#2010-08-05: JRM Mods to support variable snap block length.
#1.1 PVP Initial.


''' IMPORTS '''
import matplotlib
matplotlib.use('TkAgg')
import pylab,h5py,ratty2, time, corr, numpy, struct, sys, logging, os,ast
import iniparse
# import gps
from gps import *
import threading


# what format are the snap names and how many are there per antenna
bram_out_prefix = 'store'
# what is the bram name inside the snap block
bramName = 'bram'
verbose=True
max_level=numpy.NINF
min_level=numpy.Inf
dmax_level=numpy.NINF
dmin_level=numpy.Inf
cal_mode='full'


class GPSPoller(threading.Thread):
    def __init__(self):
        threading.Thread.__init__(self)
        self.session = gps(mode=WATCH_ENABLE)
        self.current_gps_value = None
	self._stop = threading.Event()

    def get_current_gps_value(self):
        return self.current_gps_value

    def stop(self): 
        self._stop.set()
    
    def stopped(self): 
        return self._stop.isSet()

    def run(self):
        try:
            while True:
		if self.stopped(): 
	            return
                self.current_gps_value = self.session.next()
                while not 'lat' in self.current_gps_value.keys():
                    print "Attempting to lock GPS..."
                    self.current_gps_value = self.session.next()
                    print "ok!"

        except StopIteration:
            pass


def exit_fail():
    raise
    exit()


def exit_clean():
    try:
        r.fpga.stop()
        print "Closing file."
        f.flush()
        f.close()
	gpsp.stop()
    except:
        pass
    exit()


def filewrite(stat):
    cnt=f['calibrated_spectrum'].shape[0]
    print '\nStoring entry %i...'%(cnt-1),
    sys.stdout.flush()
    for ky in stat:
        if not ky in ['raw_spectrum']: #don't store these things.
            try:
                f[ky].resize(cnt+1, axis=0)
                f[ky][cnt-1]=stat[ky]
            except KeyError:
                f.create_dataset(ky,shape=[1],maxshape=[None],data=stat[ky])
    print 'done'
    print cnt, stat['gps'], stat['gps_time'], time.ctime(stat['timestamp'])
    return cnt


def getGPS():
    try:
        # session = gps.gps("localhost", "2947")
        # session.stream(gps.WATCH_ENABLE | gps.WATCH_NEWSTYLE) 
        # report = session.next()
        # while not 'lat' in report.keys():
	#     report = session.next()
        
        report = gpsp.current_gps_value
        lat = report['lat']
	lon = report['lon']
	alt = report['alt']
	gps_time = report['time'].encode('ascii', 'ignore')

    except KeyboardInterrupt:
        exit_clean()
    # except:
    #     lat = 0.0
    #     lon = 0.0
    #     alt = 0.0
    #     gps_time = str(0.0)
    return [lat, lon, alt], str(gps_time)


def getUnpackedData(cnt):
    stat=r.get_spectrum()
    stat['gps'], stat['gps_time'] = getGPS()
    cnt=filewrite(stat)
    
    print '[%i] %s: input level: %5.2f dBm (ADC %5.2f dBm), %f degC.'%(stat['acc_cnt'],time.ctime(stat['timestamp']),stat['input_level'],stat['adc_level'],stat['adc_temp']),

    if stat['adc_shutdown']: print 'ADC selfprotect due to overrange!',
    elif stat['adc_overrange']: print 'ADC is clipping!',
    elif stat['fft_overrange']: print 'FFT is overflowing!',
    else: print 'all ok.',
    print ''
    stat['file_cnt']=cnt
    return stat




def parseargs(args):
    ret={}
    for arg in args:
        arg=arg.split('=')
        try:
            ret[arg[0]]=ast.literal_eval(arg[1])
        except ValueError:
            ret[arg[0]]=arg[1]
        except SyntaxError:
            ret[arg[0]]=arg[1]
    return ret


if __name__ == '__main__':
    from optparse import OptionParser
    p = OptionParser()
    p.set_usage('%prog [options]')
    p.add_option('-v', '--verbose', dest = 'verbose', action = 'store_true', help = 'Enable debug logging mode.')
#     p.add_option('-b', '--baseline', dest = 'baseline', action = 'store_true', default=False,
#         help = 'Keep the first trace displayed as a baseline.')
#     p.add_option('-d', '--diff', dest = 'diff', action = 'store_true', default=False,
#         help = 'Also plot the difference between the first trace and subsequent spectrum.')
#     p.add_option('-s', '--n_top', dest='n_top', type='int',default=5,
#         help='Find the top N spiky RFI candidates. Default: 5')
#     p.add_option('-f', '--play_file', dest = 'play_file', type='string', default=None,
#         help = 'Open an existing file for analysis.')
#     p.add_option('-u', '--update', dest = 'update', action = 'store_false',default=True,
#         help = 'Do not update the plots (only plot a single capture).')
    p.add_option('-p', '--skip_prog', dest='fpga_prog', action='store_false',default=True,
        help='Skip reprogramming the FPGA.')
#     p.add_option('-n', '--n_chans', dest='n_chans', type='int',default=512,
#         help='Plot this number of channels. Default: 512')
#     p.add_option('-l', '--plot_lin', dest='plot_lin', action='store_true',
#         help='Plot on linear axes. Default: semilogy.')
    p.add_option('-n', '--no_plot', dest='plot', action='store_false',default=True,
        help="Don't plot anything.")
    
    p.set_description(__doc__)
    opts, args = p.parse_args(sys.argv[1:])
    kwargs=parseargs(args)

    # n_top=opts.n_top
    verbose=opts.verbose
    # plot_baseline=opts.baseline
    # plot_diff = opts.diff
    # play_filename=opts.play_file
    cnt=0

    gpsp = GPSPoller()
    gpsp.start()


try:
    # session = gps.gps("localhost", "2947")
    # session.stream(gps.WATCH_ENABLE | gps.WATCH_NEWSTYLE) 
    # report = session.next()
    # while not 'lat' in report.keys():
    #     print "Attempting to lock GPS..."
    #     report = session.next()
    #     print "ok!"

    r = ratty2.cam.spec(**kwargs)
    co=r.cal
    print 'Config file %s parsed ok!'%(r.config['config_file'])
    print 'Connecting to ROACH %s...'%r.config['roach_ip_str'],
    r.connect()

    if verbose:
        r.logger.setLevel(logging.DEBUG)
    else:
        r.logger.setLevel(logging.INFO)
    print 'done.'

    r.initialise(skip_program=(not opts.fpga_prog), print_progress=True)
    #r.rf_frontend.stop() #disconnect from the RF interface, in case other instances want to take control while we're running.

    usrlog=('Starting file at %s (%i).'%(time.ctime(),int(time.time())))
    filename="%i.spec.h5"%(int(time.time())) 
    f = h5py.File(filename, mode="w")
    f['/'].attrs['usrlog']=usrlog

    f.create_dataset('calibrated_spectrum',shape=[1,r.config['n_chans']],maxshape=[None,r.config['n_chans']])
    f.create_dataset('gps', shape=[1,3], maxshape=[None,None])
    f.create_dataset('gps_time', data="", shape=[1,1], maxshape=[None,None])
    for key in r.config.config.keys():
        print 'Storing',key
        try:
            f['/'].attrs[key]=r.config[key]
        except:
            try:
                f[key]=r.config[key]
            except TypeError:
                if r.config[key]==None: f['/'].attrs[key]='none'
                elif type(r.config[key])==dict: 
                    f[key]=r.config[key].items()
	print f.keys()

    units='dBm'
    
    print "\n\n",usrlog,"\n\n"
    n_chans=co.config['n_chans']
    freqs=co.config['freqs']

    chan_low =co.freq_to_chan(co.config['ignore_low_freq'])
    chan_high=co.freq_to_chan(co.config['ignore_high_freq'])
    print 'Working with channels %i (%5.1fMHz) to %i (%5.1fMHz).'%(chan_low,freqs[chan_low]/1.e6,chan_high,freqs[chan_high]/1.e6)
    
    while(1):
        try:
            calData=getUnpackedData(cnt)['calibrated_spectrum']
        except KeyboardInterrupt:
	    exit_clean()

except KeyboardInterrupt:
    exit_clean()
except Exception as e:
    print e
    exit_fail()

print 'Done with all.'
exit_clean()

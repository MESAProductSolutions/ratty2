#!/usr/bin/env python

'''
Writes spectrum and GPS data from an RFI monitoring spectrometer.\n

'''
#Revisions:\n
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


import matplotlib
matplotlib.use('TkAgg')
import pylab,h5py,ratty2, time, corr, numpy, struct, sys, logging, os,ast
import iniparse
import gps

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

def exit_fail():
    raise
    exit()

def exit_clean():
    try:
        r.fpga.stop()
        print "Closing file."
        f.flush()
        f.close()
    except:
        pass
    exit()

def filewrite(stat):
    cnt=f['calibrated_spectrum'].shape[0]
    print '\nStoring entry %i...'%(cnt-1),
    sys.stdout.flush()
    #return #skip for now!
    for ky in stat:
        if not ky in ['raw_spectrum']: #don't store these things.
            try:
                #print 'Trying to store %s: '%stat,status[stat]
                f[ky].resize(cnt+1, axis=0)
                f[ky][cnt-1]=stat[ky]
            except KeyError:
                #print 'Creating dataset to store %s: '%ky,stat[ky]
                f.create_dataset(ky,shape=[1],maxshape=[None],data=stat[ky])
    print 'done'
    print cnt, stat['gps'], stat['gps_time']
    return cnt

def getGPS():
    try:
        session = gps.gps("localhost", "2947")
        session.stream(gps.WATCH_ENABLE | gps.WATCH_NEWSTYLE) 
        report = session.next()
        while not 'lat' in report.keys():
            # print "Waiting..."
	    report = session.next()
        lat = report['lat']
	lon = report['lon']
	alt = report['alt']
	gps_time = report['time'].encode('ascii', 'ignore')
    except KeyboardInterrupt:
        exit_clean()
    except:
        lat = 0.0
        lon = 0.0
	alt = 0.0
	gps_time = str(0.0)
    return [lat, lon, alt], str(gps_time)


def getUnpackedData(cnt):
    if play_filename==None:
        stat=r.get_spectrum()
	stat['gps'], stat['gps_time'] = getGPS()
        cnt=filewrite(stat)
    else:
        if cnt+1>=f['calibrated_spectrum'].shape[0]: 
            print 'No more data; end of file... bye!'
            exit_clean()
        stat={  'calibrated_spectrum':f['calibrated_spectrum'][cnt],
                'adc_overrange':f['adc_overrange'][cnt],
                'adc_shutdown':f['adc_shutdown'][cnt],
                'fft_overrange':f['fft_overrange'][cnt],
                'input_level':f['input_level'][cnt],
                'adc_level':f['adc_level'][cnt],
                'timestamp':f['timestamp'][cnt],
                'acc_cnt':f['acc_cnt'][cnt],
                'adc_temp':f['adc_temp'][cnt],
                'ambient_temp':f['ambient_temp'][cnt],
		'gps':f['gps'][cnt],
		'gps_time': f['gps_time']
                }
        cnt+=1
        #print stat['calibrated_spectrum']

    print '[%i] %s: input level: %5.2f dBm (ADC %5.2f dBm), %f degC.'%(stat['acc_cnt'],time.ctime(stat['timestamp']),stat['input_level'],stat['adc_level'],stat['adc_temp']),
    if stat['adc_shutdown']: print 'ADC selfprotect due to overrange!',
    elif stat['adc_overrange']: print 'ADC is clipping!',
    elif stat['fft_overrange']: print 'FFT is overflowing!',
    else: print 'all ok.',
    print ''
    stat['file_cnt']=cnt
    return stat

def find_n_max(data,n_max,ignore_adjacents=False):
    max_levs=numpy.ones(n_max)*-numpy.Inf
    max_locs=numpy.ones(n_max)*-numpy.Inf
    this_max_lev=-numpy.Inf
    for n,d in enumerate(data[:-1]):
        if d>this_max_lev and (not ignore_adjacents or (data[n+1]<d and data[n-1]<d)):
            #print 'Peak found at %d: %i. max now at %i'%(n,d,this_max_lev)
            loc=numpy.argmin(max_levs)
            max_levs[loc]=d
            max_locs[loc]=n
            if numpy.min(max_levs)>this_max_lev: this_max_lev=numpy.min(max_levs)
    inds = max_levs.argsort()[::-1]
    return max_levs[inds],max_locs[inds]  


# callback function to draw the data for all the required polarisations
def drawDataCallback(cnt):
    stat = getUnpackedData(cnt, None, None)
    calData=stat['calibrated_spectrum']
    
    cnt=stat['file_cnt']
    subplot1.cla()
    if stat['fft_overrange'] or stat['adc_shutdown'] or stat['adc_overrange']:
        subplot1.set_title('Spectrum %i as at %s (ADC level %5.1fdBm)'%(stat['acc_cnt'],time.ctime(stat['timestamp']),stat['adc_level']),bbox=dict(facecolor='red', alpha=0.5))
    else:
        subplot1.set_title('Spectrum %i as at %s (ADC level %5.1fdBm)'%(stat['acc_cnt'],time.ctime(stat['timestamp']),stat['adc_level']))
    subplot1.set_xlabel('Frequency (MHz)')
    subplot1.set_ylabel('Level (%s)'%units)

    if plot_baseline or plot_diff:
        subplot1.hold(True)
        subplot1.plot(freqs[chan_low:chan_high]/1.e6,baseline[chan_low:chan_high],'r',linewidth=5,alpha=0.5)

    subplot1.plot(freqs[chan_low:chan_high]/1.e6,calData[chan_low:chan_high],'b')
    #subplot1.plot(freqs[chan_low:chan_high]/1.e6,10*numpy.log10(unpackedData[chan_low:chan_high]),'b')
    #subplot1.plot(freqs[chan_low:chan_high]/1.e6,unpackedData[chan_low:chan_high],'b')

    ##collapse data for plotting:
    #collapse_factor=len(unpackedData)/plot_chans
    #collapseddata=unpackedData.reshape(plot_chans,collapse_factor).sum(1)/collapse_factor
    #if plot_type == 'lin':
    #   subplot.plot(r.freqs[::collapse_factor],collapseddata)
    #   #Plot a horizontal line representing the average noise floor:
    #   subplot.hlines((median_lev),0,r.freqs[-1])
    #elif plot_type == 'log':
    #   #subplot.semilogy(r.freqs[::collapse_factor],collapseddata)
    #   subplot.plot(r.freqs[::collapse_factor],10*numpy.log10(collapseddata))
    #   median_lev_db=10*numpy.log10(median_lev)
    #   #Plot a horizontal line representing the average noise floor:
    #   subplot.hlines(median_lev_db,0,r.freqs[-1])
    #   subplot.annotate('%3.1fdB'%(median_lev_db),(r.freqs[-1],median_lev_db))

    if plot_diff:
        dd=calData[chan_low:chan_high]-baseline[chan_low:chan_high]
        subplot2.cla()
        subplot2.plot(freqs[chan_low:chan_high]/1.e6,dd)
        subplot2.set_ylabel('Difference (dB)')
        maxs,locs=find_n_max(dd,n_top,ignore_adjacents=True)
        maxfreqs=[freqs[locs[i]+chan_low]/1.e6 for i in range(n_top)]
        for i in range(n_top):
            subplot2.annotate('%iMHz:%3.1fdB'%(numpy.round(maxfreqs[i]),maxs[i]),(maxfreqs[i],maxs[i]))
        global dmin_level
        global dmax_level
        dmin_level=min(min(dd),dmin_level)
        dmax_level=max(max(dd),dmax_level)
        subplot2.set_ylim(dmin_level-10,dmax_level+10)

    median_lev_db=numpy.median(calData[chan_low:chan_high])
    #Plot a horizontal line representing the average noise floor:
    subplot1.hlines(median_lev_db,freqs[chan_low+1]/1e6,freqs[chan_high-1]/1.e6)
    subplot1.annotate('%3.1f%s'%(median_lev_db,units),(freqs[chan_high]/1.e6,median_lev_db))

   
    #annotate:
    maxs,locs=find_n_max(calData[chan_low:chan_high],n_top,ignore_adjacents=True)
    if max(locs)==min(locs): locs=[0 for i in range(n_top)] # in case we don't find a max, locs will be [-inf, -inf, -inf...]
    maxfreqs=[freqs[locs[i]+chan_low]/1.e6 for i in range(n_top)]
    for i in range(n_top):
        print '  Local max at chan %5i (%6.2fMHz): %6.2f%s'%(locs[i]+chan_low,maxfreqs[i],maxs[i],units)
        subplot1.annotate('%iMHz:%3.1f%s'%(numpy.round(maxfreqs[i]),maxs[i],units),(maxfreqs[i],maxs[i]))

        #if plot_type == 'lin':
        #    subplot.annotate('%iMHz:%3.1fdB'%(freq,lev),(freq,collapseddata[locs[i]/collapse_factor]))
        #elif plot_type == 'log':
        #    subplot.annotate('%iMHz:%3.1fdB'%(freq,lev),(freq,10*numpy.log10(collapseddata[locs[i]/collapse_factor])))
    
    global min_level
    global max_level
    #local_min=min(calData)
    min_level=min(min(calData[chan_low:chan_high]),min_level)
    max_level=max(max(calData[chan_low:chan_high]),max_level)
    subplot1.set_ylim(min_level-10,max_level+10)
    
    fig.canvas.draw()
    if (cnt+1)<stat['acc_cnt']:
        print "WARNING: you've lost some data! Your computer's too slow for this dumprate."
    if play_filename!=None:
        print '\nPress enter to grab dataset number %i...'%cnt,
        raw_input()

    if opts.update:
        fig.canvas.manager.window.after(100, drawDataCallback, cnt)

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
    p.add_option('-b', '--baseline', dest = 'baseline', action = 'store_true', default=False,
        help = 'Keep the first trace displayed as a baseline.')
    p.add_option('-d', '--diff', dest = 'diff', action = 'store_true', default=False,
        help = 'Also plot the difference between the first trace and subsequent spectrum.')
    p.add_option('-s', '--n_top', dest='n_top', type='int',default=5,
        help='Find the top N spiky RFI candidates. Default: 5')
    p.add_option('-f', '--play_file', dest = 'play_file', type='string', default=None,
        help = 'Open an existing file for analysis.')
    p.add_option('-u', '--update', dest = 'update', action = 'store_false',default=True,
        help = 'Do not update the plots (only plot a single capture).')
    p.add_option('-p', '--skip_prog', dest='fpga_prog', action='store_false',default=True,
        help='Skip reprogramming the FPGA.')
#    p.add_option('-n', '--n_chans', dest='n_chans', type='int',default=512,
#        help='Plot this number of channels. Default: 512')
#    p.add_option('-l', '--plot_lin', dest='plot_lin', action='store_true',
#        help='Plot on linear axes. Default: semilogy.')
    p.add_option('-n', '--no_plot', dest='plot', action='store_false',default=True,
        help="Don't plot anything.")
    p.set_description(__doc__)
    opts, args = p.parse_args(sys.argv[1:])
    kwargs=parseargs(args)

    n_top=opts.n_top
    verbose=opts.verbose
    plot_baseline=opts.baseline
    plot_diff = opts.diff
    play_filename=opts.play_file
    cnt=0

try:
    if play_filename==None:
        session = gps.gps("localhost", "2947")
        session.stream(gps.WATCH_ENABLE | gps.WATCH_NEWSTYLE) 
        report = session.next()
        while not 'lat' in report.keys():
            print "Attempting to lock GPS..."
	    report = session.next()
	print "ok!"

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
            #print 'Storing',key
            try:
                f['/'].attrs[key]=r.config[key]
            except:
                try:
                    f[key]=r.config[key]
                except TypeError:
                    if r.config[key]==None: f['/'].attrs[key]='none'
                    elif type(r.config[key])==dict: 
                        f[key]=r.config[key].items()
#                        print 'Stored a dict!'
	print f.keys()
    else:
        print 'Opening file %s...'%play_filename
        f=h5py.File(play_filename,'r')
        conf_ovr=dict(f['/'].attrs)
        for key in f.keys():
            if not key in ['calibrated_spectrum','timestamp','adc_overrange','fft_overrange','adc_shutdown','adc_level','input_level']:
                #print 'trying',key
                if len(f[key])>1: conf_ovr[key]=f[key][:]
                else: conf_ovr[key]=f[key]

        co=ratty2.cal.cal(**conf_ovr)
        usrlog=f['/'].attrs['usrlog']

    if co.config['antenna_bandpass_calfile'] == 'none':
        units='dBm'
    else:
        units='dBuV/m'

    print "\n\n",usrlog,"\n\n"
    n_chans=co.config['n_chans']
    freqs=co.config['freqs']

    bp=co.config['system_bandpass']
    af=co.config['ant_factor']

    chan_low =co.freq_to_chan(co.config['ignore_low_freq'])
    chan_high=co.freq_to_chan(co.config['ignore_high_freq'])
    print 'Working with channels %i (%5.1fMHz) to %i (%5.1fMHz).'%(chan_low,freqs[chan_low]/1.e6,chan_high,freqs[chan_high]/1.e6)
    if play_filename != None:
        # set up the figure with a subplot for each polarisation to be plotted
        fig = matplotlib.pyplot.figure()
        if opts.diff or opts.baseline:
            baseline=numpy.zeros(n_chans)
            print 'Fetching baseline...',
            sys.stdout.flush()
            baseline=getUnpackedData(cnt, None, None)['calibrated_spectrum']
            print 'done'
        if opts.diff: 
            subplot1 = fig.add_subplot(2, 1, 1)
            subplot2 = fig.add_subplot(2, 1, 2)
        else: subplot1 = fig.add_subplot(1, 1, 1)
        fig.canvas.manager.window.after(100, drawDataCallback,cnt)
        matplotlib.pyplot.show()
        print 'Plot started.'
    else:
    	while(1):
            try:
                calData=getUnpackedData(cnt)['calibrated_spectrum']
                # maxs,locs=find_n_max(calData[chan_low:chan_high],n_top,ignore_adjacents=True)
                # maxfreqs=[freqs[locs[i]+chan_low]/1.e6 for i in range(n_top)]
                # for i in range(n_top):
                #     print '  Local max at chan %5i (%6.2fMHz): %6.2f%s'%(locs[i]+chan_low,maxfreqs[i],maxs[i],units)
	    except KeyboardInterrupt:
	        exit_clean()

except KeyboardInterrupt:
    exit_clean()
except Exception as e:
    print e
    exit_fail()

print 'Done with all.'
exit_clean()

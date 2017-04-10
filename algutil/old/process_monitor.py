#!/usr/bin/python


import subprocess
import sys
import resource
import collections
import datetime
import socket


def execute_str(cmd_str):

    p = subprocess.Popen(cmd_str,shell=True,stdout = subprocess.PIPE,stderr = subprocess.PIPE)
    (stdout,stderr) = p.communicate()
    return_code = p.returncode
    err = return_code!=0
    stdout_str = str(stdout)
    stderr_str = str(stderr)

    return (err,stdout_str,stderr_str)

def write_msg_to_file(filepath,msg):
    fid = open(filepath,'wt')
    fid.write(msg+'\n')
    fid.close()


def get_timestamp():
    t=datetime.datetime.now()
    try:
        #Python 3 version
        timestamp='{0.year:d}_{0.month:02d}_{0.day:02d}__{0.hour:02d}_{0.minute:02d}_{0.second:02d}'.format(t)
    except:
        #Python 2 version
        timestamp = '%d_%02d_%02d__%02d_%02d_%02d'%(t.year,t.month,t.day,t.hour,t.minute,t.second)
    return timestamp


def get_timestamp_delta(ts_begin,ts_end):
    begin_secs = timestamp_to_seconds(ts_begin)
    end_secs = timestamp_to_seconds(ts_end)
    duration = end_secs - begin_secs
    return duration

def timestamp_to_seconds(timestamp):
    (year,month,day,junk,hour,minute,second)= timestamp.split('_')
    (year,month,day,junk,hour,minute,second) = (
        int(year),int(month),int(day),junk,int(hour),int(minute),int(second))

    dt = datetime.datetime(year,month,day,hour,minute,second)
    days = dt.toordinal()
    days = days - 733000 # keep the seconds to a managable size...
    secs = days*24*60*60 + hour*60*60 + minute*60 + second
    return secs




cmd_str = ' '.join(sys.argv[1:])




(err,stdout_str1,stderr_str) = execute_str('whoami')
(err,stdout_str2,stderr_str) = execute_str('uname -a')
(err,stdout_str3,stderr_str) = execute_str('df -h')
(err,stdout_str4,stderr_str) = execute_str('echo num processors: `grep -c ^processor /proc/cpuinfo`')
(err,stdout_str5,stderr_str) = execute_str('free -th')

(err,stdout_str20,stderr_str) = execute_str('pwd')


(err,stdout_str10,stderr_str) = execute_str('curl "http://metadata.google.internal/computeMetadata/v1/project/project-id" -H "Metadata-Flavor: Google"')
(err,stdout_str11,stderr_str) = execute_str('curl "http://metadata.google.internal/computeMetadata/v1/instance/zone" -H "Metadata-Flavor: Google"')
(err,stdout_str12,stderr_str) = execute_str('curl "http://metadata.google.internal/computeMetadata/v1/instance/machine-type" -H "Metadata-Flavor: Google"')
(err,stdout_str13, stderr_str) = execute_str('curl "http://metadata.google.internal/computeMetadata/v1/instance/hostname" -H "Metadata-Flavor: Google"')
(err,stdout_str14, stderr_str) = execute_str('curl "http://metadata.google.internal/computeMetadata/v1/instance/scheduling/preemptible" -H "Metadata-Flavor: Google"')
(err,stdout_str15, stderr_str) = execute_str('curl "http://metadata.google.internal/computeMetadata/v1/instance/description" -H "Metadata-Flavor: Google"')
# todo - add info about gce disk types and sizes via instance/disks

(err,stdout_str16, stderr_str) = execute_str('curl "http://metadata.google.internal/computeMetadata/v1/instance/disks/" -H "Metadata-Flavor: Google" | wc -l')
num_disks = int(stdout_str16)
disk_info = ''
for d in range(num_disks):
    (err, stdout_str17, stderr_str) = execute_str(
        'curl "http://metadata.google.internal/computeMetadata/v1/instance/disks/%d/type" -H "Metadata-Flavor: Google" | wc -l')
    info = 'disk %d: %s; '%(d,stdout_str17)
    disk_info += info

#cmd_str

#write to log
fid = open('monitor_status_before.txt','w')
fid.write(stdout_str1)
fid.write(stdout_str2)
fid.write(stdout_str3)
fid.write(stdout_str4)
fid.write(stdout_str5)


fid.write('project-id: %s'%stdout_str10)
fid.write('zone: %s'%stdout_str11)
fid.write('machine-type: %s'%stdout_str12)
fid.write('hostname: %s'%stdout_str13)
fid.write('preemptable: %s'%stdout_str14)
fid.write('description: %s'%stdout_str15)
fid.write('disk info: %s'%disk_info)


fid.write('\n\ncd %s\n\n%s\n\n'%(stdout_str20, cmd_str))
fid.close()


#write to stdout
print(stdout_str1)
print(stdout_str2)
print(stdout_str3)
print(stdout_str4)
print(stdout_str5)


print('project-id: %s'%stdout_str10)
print('zone: %s'%stdout_str11)
print('machine-type: %s'%stdout_str12)
print('hostname: %s'%stdout_str13)
print('preemptable: %s'%stdout_str14)
print('description: %s'%stdout_str15)
print('disk info: %s'%disk_info)

print('\n\ncd %s\n\n%s\n\n' % (stdout_str20, cmd_str))



# launch dstat as a background process
dstat_cmd_str = "dstat --nocolor --noheaders --freespace --top-cpu-adv --top-io --top-mem --top-bio-adv --dstat-mem -cdngymt --output dstat.log.txt"
dstat_list=dstat_cmd_str.split()
dstat_process = subprocess.Popen(dstat_list)
#df_monitor_cmd_str="/home/disk_monitor.sh df.log.txt"
#df_monitor_process=subprocess.Popen(df_monitor_cmd_str.split())



start_timestamp = get_timestamp()



# run command, wait it to complete
exit_code = subprocess.call(cmd_str, shell=True)

end_timestamp = get_timestamp()

# kill the dstat background process in a way that allows it to close the file properly
dstat_process.terminate()
#df_monitor_process.terminate()

# generate stats from running
usage = resource.getrusage(resource.RUSAGE_CHILDREN)
if exit_code == 0:
    passfail = "pass"
else:
    passfail = "FAIL"
maxrss_memory = float(usage.ru_maxrss) / (2 ** 20)  # convert memory to GB

usage_dict = collections.OrderedDict()

usage_dict['#status'] = passfail
usage_dict['exit_code'] = str(exit_code)
usage_dict['start_time'] = start_timestamp
usage_dict['end_time'] = end_timestamp
usage_dict['wallclock_duration_s'] = str(get_timestamp_delta(start_timestamp, end_timestamp))
usage_dict['user_cputime_s'] = '%7.3f' % usage.ru_utime
usage_dict['system_cputime_s'] = '%7.3f' % usage.ru_stime
# note maxrss may just capture the real memory used by the largest child, not all the children
usage_dict['maxrss_memory_gb'] = '%7.3f' % maxrss_memory
usage_dict['page_faults'] = '%d' % usage.ru_majflt
usage_dict['node_ip_addr'] = socket.gethostbyname(socket.gethostname())


(err,stdout_str3,stderr_str) = execute_str('df -h')
(err,stdout_str5,stderr_str) = execute_str('free -th')


#write to log
fid = open('monitor_status_after.txt','a')

fid.write(stdout_str3)
fid.write(stdout_str5)

for key in usage_dict:
    fid.write('%s: %s\n'%(key,usage_dict[key]))
fid.close()


# write to stdout
print(stdout_str3)
print(stdout_str5)

for key in usage_dict:
    print('%s: %s' % (key, usage_dict[key]))

sys.exit(exit_code)


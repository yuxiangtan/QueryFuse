#   Copyright {2015} Yuxiang Tan
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import os, sys, shutil, time, shlex, subprocess, random, getopt
#import json

def checkParameter(param,key,dType,allowed=[],checkFile=False,optional=False):
    #generic function that checks if a parameter was in the parameter file, casts to the right data type and if the parameter is a file/ directory checks if it actually exists
    if optional and key not in param:
        param[key]=''
    else:
        #check if key is present
        if key not in param: print 'Parameter '+key+' is missing in the parameter file.'; sys.exit(1)

        #cast to correct data type
        if dType == bool:
            param[key] = param[key] in ['True','TRUE','true','T','1']
        else:
            param[key]=dType(param[key])

        #check if the value for the key is allowed
        if len(allowed)>0 and param[key] not in allowed:
            print 'Parameter '+key+' can only be one of the following: '+allowed; sys.exit(1)
    
        #if file or directory check if it exists
        if checkFile and not os.path.exists(param[key]):
            print 'The file in '+key+' = ',param[key],' does not exist'; sys.exit(1)
            
#This function constructs a dictionary with pairs: parameter name - value
def parse_parameters(par_file):
    dict_param = dict([(l.split("=")[0].strip(), l.split("#")[0].split("=")[1].strip()) for l in par_file.readlines() if "=" in l and l.strip()[0]!="#"])
    return dict_param

def update_parameters(args):
    #This function updates the parameter dictionaries with additional values provided in the command line
    pnames = args[0::2]
    tuple_args = zip(pnames, args[1::2])
    file_param = filter(lambda x: x[0]== '-p', tuple_args)
    # read parameter file
    fname = file_param[0][1]
    param = parse_parameters(open(fname))
    # read additional parameters
    other_param = filter(lambda x: x[0]!= '-p', tuple_args)
    new_param = []
    if len(other_param) > 0:
        # check additional parameter names
        for t_arg in other_param:
            if '--' in t_arg[0] and t_arg[0].split('--')[1] in param.keys():
                new_param.append((t_arg[0].split('--')[1], t_arg[1]))
            else:
                print "Parameter %s is specified incorrectly, please see parameter file for all possible parameters." %t_arg[0]
                sys.exit(0) 
    # update parameters
    new_dict = dict([t for t in new_param if t[0] in param.keys()])
    param.update(new_dict)
    return param, fname, new_dict

#this function need to be personalized into this pipeline!
def initialize_standard(param):
    #check default pipeline parameters and create working directories
    checkParameter(param,key='raw_filenames',dType=str,checkFile=True)
    checkParameter(param,key='accepted_bam',dType=str,checkFile=True)
    checkParameter(param,key='unmapped_bam',dType=str,checkFile=True)
    checkParameter(param,key='paired',dType=bool)
    checkParameter(param,key='clean_run',dType=bool)
    checkParameter(param,key='verbose',dType=bool)
    checkParameter(param,key='raw_file_header',dType=bool)
    checkParameter(param,key='whole_gene_list',dType=str,checkFile=True)
    checkParameter(param,key='genome_fa',dType=str,checkFile=True)
    checkParameter(param,key='tophat_genome_ref',dType=str,checkFile=False)
    checkParameter(param,key='read_len',dType=int)
    checkParameter(param,key='read_std',dType=int)
    checkParameter(param,key='Align_percent',dType=int)
    checkParameter(param,key='split_n',dType=int)
    checkParameter(param,key='span_n',dType=int)
    checkParameter(param,key='sum_n',dType=int)
    checkParameter(param,key='num_processors',dType=int)
    checkParameter(param,key='timeout',dType=int)
    checkParameter(param,key='span_only_filter',dType=int)
    checkParameter(param,key='step_size_query',dType=str)
    checkParameter(param,key='step_size_other',dType=str)
    
    #checking working directory and going there
    checkParameter(param,key='working_dir',dType=str,checkFile=False) 
    if  param['working_dir'][len(param['working_dir'])-1]!='/':
        param['working_dir']=param['working_dir']+'/'
    #  os.chdir(os.path.dirname(param['working_dir']))

    #directory where all the scripts are located
    checkParameter(param,key='scripts_dir',dType=str)
    if  param['scripts_dir'][len(param['scripts_dir'])-1]!='/':
        param['scripts_dir']=param['scripts_dir']+'/'
    checkParameter(param,key='tophat_genome_ref',dType=str)
    if not os.path.exists(param["tophat_genome_ref"]+".1.bt2"):
	print 'tophat_genome_ref does not exist'; sys.exit(0) 
    param['no_module_load']="TRUE"
    if param['num_processors']==1:
        param['run_single_cpu']="TRUE"
    else:
	param['run_single_cpu']="FALSE"
    
    #if directory exists and the pipeline should be run from scratch delete the directory
    if param['clean_run']:
        #check before deleting any previous results
        answer = input('Are you sure you want to delete all existing results? (yes/no): ')
        if answer != 'yes':
            print 'Stopping... If you want to resume, please adjust the parameter file.'
            exit(0)
        if os.path.exists(param['working_dir']+'results/'):
            shutil.rmtree(param['working_dir']+'results/')
        #not doing this yet, because the resume function is still not working well.
	#if os.path.exists(param['working_dir']+'bams/'):
        #    shutil.rmtree(param['working_dir']+'bams/')
        if os.path.exists(param['working_dir']+'intermedias/'):
            shutil.rmtree(param['working_dir']+'intermedias/')

    #if results or report directory do not exist create them
    if not os.path.exists(param['working_dir']+'results/'):
         os.makedirs(param['working_dir']+'results/')
    if not os.path.exists(param['working_dir']+'bams/'):
         os.makedirs(param['working_dir']+'bams/')
    if not os.path.exists(param['working_dir']+'intermedias/'):
         os.makedirs(param['working_dir']+'intermedias/')
	 
    #hard code min_score here
    param["min_score"]="11"

def write_updated_file(updated, param, parameter_file):
    new_parameter_file=param['working_dir']+'results/'+parameter_file[:-4]+'_used.txt'
    with open(parameter_file) as f:
        with open(new_parameter_file, "w") as f1:
            for l in f:
                k = l.split("=")[0].strip()
                if "=" in l and k in updated.keys():
                    f1.write('%s= %s\n' %(l.split("=")[0], updated[k]))
                else:
                    f1.write(l)
    param['parameter_file']=new_parameter_file
   
def readQueryIDFilenames(param):
#Get filenames and stubs, while checking if the files actually exist in paired mode there will be 2 fastq files per sample, otherwise only 1
    param['stub']=[]
    param['raw_files']=[]
    #if param['paired']: param['raw_files2']=[]
    
    #assign raw file locations
    f = open (param['raw_filenames'] , 'r')
    #skip first line if there is a header specified
    header=param['raw_file_header']
    for line in f:
        if header: 
            header=False
        else:
            line = line.strip()
            line = line.split('\t')
            #param['raw_files'].append(line[0])
            ## check if filename actually exists
            #if not os.path.exists(line[0]):
            #    print 'The file '+line[0]+' does not exist';
            #    sys.exit(0)
            param['stub'].append(line[0])
            
    param['num_samples']=len(param['stub'])
    
    #assign the raw files the fastq_file list - this is a convenience ideally you should speficfy that the pipeline starts with teh raw_files
    #param['QueryID_files']=param['raw_files'][:]
       
    #start a log file that keeps track on which files are successfully completed
    param['run_log']=[[True]*param['num_samples']]
    param['run_log_headers']=['raw']    

def initialize_logfiles(param):
    ##create log files and log directory if they do not exist already
    #log_dir=param['working_dir']+'results/log/'
    ##if the doesn't exist create it
    #if not os.path.exists(log_dir):
    #     os.makedirs(log_dir)
    ##if log files don't exit create them
    #for idx in range(param['num_samples']):
    #    log_error_file = log_dir+param['stub'][idx]+'_error.log'
    #    log_out_file = log_dir+param['stub'][idx]+'_out.log'
    #    if not os.path.exists(log_error_file):
    #        open(log_error_file,'a').close()
    #    if not os.path.exists(log_out_file):
    #        open(log_out_file,'a').close()    
    param['log_handle'] = param['working_dir']+'results/main.log'

#writeLog will just write in the main log, but not in others sub logs.
def writeLog(string,param):
    if param['verbose']: print string
    handle=open(param['log_handle'],'a')
    handle.write(string)
    handle.close()

def addInputOutputFiles(param,input_files,output_files):
    #store information on which file type to work and what output file to use
    param['input_files'] = input_files
    param['output_files'] = output_files

    #check if the specified input files actually exist otherwise throw an error
    if not param.has_key(param['input_files']):
        writeLog('The input file type you specified does not exist')
        sys.exit(0)
        
    #check if the key for the specified output files actually exist otherwise create them
    if not param.has_key(param['output_files']) and param['output_files']!='':
         param[param['output_files']]=['']*param['num_samples']

def check_if_queue_finished_sucessful(param):
  # check how many ended successfully
    s_noend = []
    for idx in range(param['num_samples']):
        logfile = open(param['working_dir']+'results/log/'+param['stub'][idx]+'_out.log')
        lines_end = [line for line in logfile.readlines() if 'ENDING %s |' %(param['current_flag']) in line.rstrip()]
        #check if the job was successfully finished
        if len(lines_end) == 1:
            #check if there is suppose to be a output file and if that is the case store the location of the output file
            if  param['output_files']!='':
                retval=lines_end[0].split('|')[1].strip()
                if len(retval.split(';'))==1:
                    param[param['output_files']][idx]=retval
                else:
                    param[param['output_files']][idx]=retval.split(';')[0]
                    param[param['output_files']+'2'][idx]=retval.split(';')[1]    
            #add an entry into the run log that shows that the current job has been finished successfully
            param['run_log'][-1][idx]=True
        else:
            if  param['output_files']!='': param[param['output_files']][idx]=''
            s_noend.append(param['stub'][idx])
            
    sucessful=True
    #output error if not all samples finished successfully otherwise indicate that this step was completed sucessfully
    if len(s_noend) > 0:
        writeLog('error in samples %s' %(';'.join([s for s in s_noend])),param)
        writeLog('\n',param)
        sucessful=False
    else:        
        writeLog(param['current_flag']+' successful!\n\n',param)
    return sucessful


def dumpParameters(param):
    #dumps the parameter file into a JSON object
    with open(param['working_dir']+'results/parameters.json', 'w') as f:
        json.dump(param, f)

def resume_func(cmd, resume_stat, step_name, next_step_name, LOG_OUT):
    #may I should add log_out.error here.
    if resume_stat == 1:
        h_LOG_OUT=open(LOG_OUT,"r")
        count_res = 0
        for line in h_LOG_OUT:
	    if (step_name+" done") in line:
                count_res +=1
	h_LOG_OUT.close()            
        if count_res==0:
            #subprocess.call(cmd,shell=True)
	    status=sync_cmd(cmd)
	#    for i in status:
	#	print(i)
	    stdout_re=status[0]
	    stderr_re=status[1]
	    retcode_re=status[2]
	    resume_stat=0
        else:
	    stdout_re=""
	    stderr_re=""
	    retcode_re=0
            print step_name+" pass"
    else:
        log_outf=open(LOG_OUT,"a")
        #subprocess.call(cmd,shell=True)
	status=sync_cmd(cmd)
	#for i in status:
	#    print(status)
	stdout_re=status[0]
	stderr_re=status[1]
        retcode_re=status[2]
	if retcode_re==0:
	    log_outf.write(step_name+" done\n")
	else:
	    log_outf.write(step_name+" failed\n")
	    
        log_outf.close()
    return (resume_stat,stdout_re, stderr_re, retcode_re)

def sync_cmd(cmd):    
    proc = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    retcode = proc.returncode
    return (stdout, stderr, retcode)

#need to be optimize to my purpose
def initialize_module():
    #Import modules
    import sys, getopt, json

    working_dir='./'
    #Check arguments
    if len(sys.argv)<5:
        print 'ERROR: Specify the index of the file the parameter should be run on.'
        sys.exit(0)
    optlist,cmdlist = getopt.getopt(sys.argv[1:],'i:n:d:')
    for opt in optlist:
        if opt[0]=='-i': file_index=opt[1]
        if opt[0]=='-n': num_processors=opt[1]
        if opt[0]=='-d': working_dir=opt[1]
    print '###########################'
    print sys.argv
    print working_dir

    #Read and initialize parameters
    with open(working_dir+'results/parameters.json') as f:
        param = json.load(f)
    param['file_index']=int(file_index)
    param['num_processors']=num_processors
    
    #use the input files that were specified in the pipeline call
    #param['working_file']=param[param['input_files']][param['file_index']]
    #if param['paired'] and param['input_files']+'2' in param:
    #   param['working_file2']=param[param['input_files']+'2'][param['file_index']]

    #name of the log file but it is not really the way I want it to be.
    log_file = param['working_dir']+'results/log/'+param['stub'][param['file_index']]+"/"+param['stub'][param['file_index']]+'.log'

    return param
     
def key_step_check(cmd_status, step_name, log_whole, log_error):
    if cmd_status[3]==0:
        log_whole.write(cmd_status[1]+cmd_status[2])
    else:
        log_error.write(cmd_status[2])
        log_whole.write(cmd_status[1]+cmd_status[2])
        sys.exit('Warning: '+step_name+' failed, please see the error.log of that section for more details.')
	
def optional_step_check(cmd_status, log_whole, log_error, shared_str):
    if cmd_status[3]==0:
        log_whole.write(cmd_status[1]+cmd_status[2])
	print "finished "+shared_str
        log_whole.write("finished "+shared_str+'\n')
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n')
        log_whole.write("===================="+'\n')
    else:
        log_error.write(cmd_status[2])
        log_whole.write(cmd_status[1]+cmd_status[2])
	print "Warning: failed to "+shared_str
        log_whole.write("Warning: failed to "+shared_str+'\n')
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n')
        log_whole.write("===================="+'\n')

#read fa file into dictionary
def fa_to_dic(fa_file, row_distant):
	h_FA=open(fa_file,"r")
	fa_dic={}
	while 0==0:
		fa_line=h_FA.readline()
		if fa_line=="":
			break
		fa_line_ID=fa_line.strip().split(">")[1]
		seq_line=""
		for i in range(row_distant):
			seq_line+=h_FA.readline().strip()
		fa_dic[fa_line_ID]=seq_line
	h_FA.close()
	return(fa_dic)

#to get the row distant between two reads in a fa file, in order to help read the fa correctly no matter how many lines are in between.
def fa_to_row_distance(fa_file,log_error):
    h_FA=open(fa_file,"r")
    FA_line1=h_FA.readline()
    if FA_line1[0]!=">":
	print "Warning: The input of unmapped.bam_first_mate.fa in not in correct format with > at the first row in QF_summary_process.py, exit.\n, exit."
	log_error.write("The input of unmapped.bam_first_mate.fa in not in correct format with > at the first row in QF_summary_process.py, exit.\n"); sys.exit(1)
    row_int_fa=0
    for line in h_FA:
	if line[0]==">":
	    break
	else:
	    row_int_fa+=1
    h_FA.close()
    return(row_int_fa)

#read the (subtract) bed file into a dictionary, it is much more slower than I expected, maybe it is the function to improve speed.
def subtract_bed_to_dic(bed_file):
    h_subtract_bed=open(bed_file,"r")
    bed_dic={}
    for line in h_subtract_bed:
	line_row=line.strip().split("\t")
	#in case it has the end num annotation
	line_ID=line_row[0].split("/")[0]
	if not line_ID in bed_dic.keys():
	    bed_dic[line_ID]={}
	    #this is to eliminate duplicated lines.
	    bed_dic[line_ID]["^".join(line_row)]=""
	else:
	    bed_dic[line_ID]["^".join(line_row)]=""
    h_subtract_bed.close()
    return(bed_dic)
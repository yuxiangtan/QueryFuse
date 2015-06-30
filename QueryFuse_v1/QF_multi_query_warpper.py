# -*- coding: cp936 -*-

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

"""
Usage: python QF_multi_query_wrapper.py -p parameter_file.txt
-h help
-p parameter file				*[No default value]

============================

This script is the wrapper of QueryFuse.
It reads the parameter file and recognize and check the parameters.
Instead of using qsub system, it will seperate jobs into small modules and organize their running in a node with accessable processors.
This is the merge version of QF_wrapper and QF_each_query.
============================

Python & Module requirement:
Versions: 2.7 or above
Module: No additional Python Module is required.

============================
Library file requirement:
Not Standalone version, few library file is required.
bowtie2, samtools, bedtools are need to reinstall and preloaded by anaconda. 
============================

"""
##intial functions
def run_proc_prepare(param,query,result_outdir):
    out_log=sync_cmd(cmd_prepare_for_query(param,query,result_outdir))
    return(out_log)

def run_proc_pair(param,query):
    out_log=async_cmd(cmd_paired_end_for_query(param,query))
    return(out_log)

def run_proc_single(param,query):
    out_log=async_cmd(cmd_single_end_for_query(param,query))
    return(out_log)

def each_query(param):
    # Extract real values from param dict
    num_parallels = param['num_processors']
    print "num_parallels: "+str(num_parallels)
    queries = param['stub']
    timeout = int(param['timeout'])
        
    processes = {}
    proc_prepare={}
    results = {}
    cur_num_procs = 0
    proc_summary = {}
    # Dispatch each query
    if num_parallels>1:
        for query in queries:
            #setup parameters
            #run create output directories (mainly will go into two divisions: intermediate and result.)
            #print to test for leak
            
            BAM_dir = param['bam_dir']
            inter_dir = param['file_prefix']
            result_dir = param['result_folder_prefix']
            
            #BAM_outdir = param['bam_dir']+param['stub'][param['file_index']]+'/'
            inter_outdir = param['file_prefix']+query+'/'
            result_outdir = param['result_folder_prefix']+query+'/'
            log_outdir = param['LOG_F']+query+'/'
            
            #run_log?????? preprocess (it is not at this level), pair, single, summary each have one run log? And a log for each of these big step?
            log_outdir_prepare_for_query = log_outdir+"log_of_prepare_for_query/"
            log_outdir_paired_end = log_outdir+"log_of_paired_end_run/"
            log_outdir_single_end = log_outdir+"log_of_singleton_unmapped_run/"
            log_outdir_summary = log_outdir+"log_of_summary_run/"
            
            #if not os.path.exists(BAM_outdir):
            #    os.makedirs(BAM_outdir)
            #check path and create folders
            if not os.path.exists(inter_outdir):
                os.makedirs(inter_outdir)
            
            if not os.path.exists(result_outdir):
                os.makedirs(result_outdir)
                
            if not os.path.exists(log_outdir):
                os.makedirs(log_outdir)
            
            if not os.path.exists(log_outdir_prepare_for_query):
                os.makedirs(log_outdir_prepare_for_query)
            
            if not os.path.exists(log_outdir_paired_end):
                os.makedirs(log_outdir_paired_end)
            
            if not os.path.exists(log_outdir_single_end):
                os.makedirs(log_outdir_single_end)
               
            if not os.path.exists(log_outdir_summary):
                os.makedirs(log_outdir_summary)
            #print query
            #print time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))
            #print cur_num_procs
            #
            # If there are not enough free processes
            while int(cur_num_procs) + 2 > int(num_parallels):
                # Wait for enough processes to terminate
                #print "Hold at the waiting line"
                cur_num_procs -= wait_for(processes, results, timeout)
            # Launch jobs for this query
            proc_prepare[query] = run_proc_prepare(param,query,result_outdir)
            h_prepare_OUT=open(log_outdir_prepare_for_query+"out.log","w")
            h_prepare_OUT.write(proc_prepare[query][0]+proc_prepare[query][1])
            h_prepare_OUT.close()
            if proc_prepare[query][2]!=0:
                QF_all_modules.writeLog('prepare for '+query+' failed in each_query module\n',param)
                QF_all_modules.writeLog(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+" \n",param)
                h_prepare_ERR=open(log_outdir_prepare_for_query+"error.log","w")
                h_prepare_ERR.write(proc_prepare[query][1])
                h_prepare_ERR.close()
                #since the praparing failed, the down stream will not be submitted for this query
                continue
            else:
                QF_all_modules.writeLog('finished preparing for '+query+' in each_query module\n',param)
                QF_all_modules.writeLog(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+" \n",param)
            proc_pair = run_proc_pair(param,query)
            proc_single = run_proc_single(param,query)
            # Store jobs about this query in a dictionary
            jobs = dict(single = proc_single, pair = proc_pair)
            processes[query] = jobs
            # Increment number of processes active
            cur_num_procs += 2
        # Wait for all extant jobs to finish
        while(cur_num_procs != 0):
            #print "Wait for all the rest job to finish"
            cur_num_procs -= wait_for(processes, results, timeout)
        #the following loop is replaced by the logging function in wait_for module.
        #for queries in results.keys():
        #    print queries
        #    for jobs in results[queries].keys():
        #        print jobs
        #        for logs in results[queries][jobs]:
        #            if type(logs) is types.IntType:
        #                print logs
        #            else:
        #                print logs[0]
    else:
        #single node mode:
        proc_pair={}
        proc_single={}
        for query in queries:
            #setup parameters
            #run create output directories (mainly will go into two divisions: intermediate and result.)
            BAM_dir = param['bam_dir']
            inter_dir = param['file_prefix']
            result_dir = param['result_folder_prefix']
            
            #BAM_outdir = param['bam_dir']+param['stub'][param['file_index']]+'/'
            inter_outdir = param['file_prefix']+query+'/'
            result_outdir = param['result_folder_prefix']+query+'/'
            log_outdir = param['LOG_F']+query+'/'
            
            #run_log?????? preprocess (it is not at this level), pair, single, summary each have one run log? And a log for each of these big step?
            log_outdir_prepare_for_query = log_outdir+"log_of_prepare_for_query/"
            log_outdir_paired_end = log_outdir+"log_of_paired_end_run/"
            log_outdir_single_end = log_outdir+"log_of_singleton_unmapped_run/"
            log_outdir_summary = log_outdir+"log_of_summary_run/"
            
            #if not os.path.exists(BAM_outdir):
            #    os.makedirs(BAM_outdir)
            #check path and create folders
            if not os.path.exists(inter_outdir):
                os.makedirs(inter_outdir)
            
            if not os.path.exists(result_outdir):
                os.makedirs(result_outdir)
                
            if not os.path.exists(log_outdir):
                os.makedirs(log_outdir)
            
            if not os.path.exists(log_outdir_prepare_for_query):
                os.makedirs(log_outdir_prepare_for_query)
            
            if not os.path.exists(log_outdir_paired_end):
                os.makedirs(log_outdir_paired_end)
            
            if not os.path.exists(log_outdir_single_end):
                os.makedirs(log_outdir_single_end)
               
            if not os.path.exists(log_outdir_summary):
                os.makedirs(log_outdir_summary)
            
            #Because this has only one processer, jobs will be run one by one, and saved in different format in a dictionary
            proc_prepare[query] = sync_cmd(cmd_prepare_for_query(param,query,result_outdir))
            h_prepare_OUT=open(log_outdir_prepare_for_query+"out.log","w")
            h_prepare_OUT.write(proc_prepare[query][0]+proc_prepare[query][1])
            h_prepare_OUT.close()
            if proc_prepare[query][2]!=0:
                QF_all_modules.writeLog('prepare for '+query+' failed in each_query module\n',param)
                QF_all_modules.writeLog(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+" \n",param)
                h_prepare_ERR=open(log_outdir_prepare_for_query+"error.log","w")
                h_prepare_ERR.write(proc_prepare[query][1])
                h_prepare_ERR.close()
                #since the praparing failed, the down stream will not be submitted for this query
                continue
            else:
                QF_all_modules.writeLog('finished preparing for '+query+' in each_query module\n',param)
                QF_all_modules.writeLog(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+" \n",param)
            proc_pair[query] = sync_cmd(cmd_paired_end_for_query(param,query))
            #since the pair_end and single_end are parallel and not necessary to have input reads, now I will just report failed in the mail log.
            #h_pair_end_OUT=open(log_outdir_paired_end+"out.log","w")
            #h_pair_end_OUT.write(proc_pair[query][0]+proc_pair[query][1])
            #h_pair_end_OUT.close()
            if proc_pair[query][2]!=0:
                QF_all_modules.writeLog('pair process in :'+query+' failed in each_query module\n',param)
                QF_all_modules.writeLog(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+" \n",param)
                #h_pair_end_ERR=open(log_outdir_paired_end+"error.log","w")
                #h_pair_end_ERR.write(proc_pair[query][1])
                #h_pair_end_ERR.close()
            else:
                QF_all_modules.writeLog('finished processing pair in :'+query+' in each_query module\n',param)
                QF_all_modules.writeLog(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+" \n",param)
            proc_single[query] = sync_cmd(cmd_single_end_for_query(param,query))
            #h_singleton_OUT=open(log_outdir_single_end+"out.log","w")
            #h_singleton_OUT.write(proc_single[query][0]+proc_single[query][1])
            #h_singleton_OUT.close()
            if proc_single[query][2]!=0:
                QF_all_modules.writeLog('single process in :'+query+' failed in each_query module\n',param)
                QF_all_modules.writeLog(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+" \n",param)
                #h_singleton_ERR=open(log_outdir_single_end+"error.log","w")
                #h_singleton_ERR.write(proc_single[query][1])
                #h_singleton_ERR.close()
            else:
                QF_all_modules.writeLog('finished processing single in :'+query+' in each_query module\n',param)
                QF_all_modules.writeLog(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+" \n",param)
            proc_summary[query] = sync_cmd(cmd_summary_end_for_query(param,query))
            #h_summary_OUT=open(log_outdir_summary+"out.log","w")
            #h_summary_OUT.write(proc_summary[query][0]+proc_summary[query][1])
            #h_summary_OUT.close()
            if proc_summary[query][2]!=0:
                QF_all_modules.writeLog('summary process in :'+query+' failed in each_query module\n',param)
                QF_all_modules.writeLog(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+" \n",param)
                #h_summary_ERR=open(log_outdir_summary+"error.log","w")
                #h_summary_ERR.write(proc_summary[query][1])
                #h_summary_ERR.close()
            else:
                QF_all_modules.writeLog('finished processing summary in :'+query+' in each_query module\n',param)
                QF_all_modules.writeLog(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+" \n",param)
            h_overall_log=open(param['LOG_F']+query+'/'+"whole_procedure.log","a")
            h_overall_log.write('pair\n'+proc_pair[query][0]+'\n'+proc_pair[query][1]+'\n')
            h_overall_log.write('single\n'+proc_single[query][0]+'\n'+proc_pair[query][1]+'\n')
            h_overall_log.write('summary\n'+proc_summary[query][0]+'\n'+proc_pair[query][1]+'\n')
            h_overall_log.write
    #print(processes)
    #print(results)
    

def cmd_prepare_for_query(param,query,result_outdir):
    GENE_ID=query
    GENE_BED=param['whole_gene_list']
    QUERY_BED=result_outdir+"query_gene.bed"
    QUERY_FA=result_outdir+"query_gene.fa"
    cmd="echo 'get the bed file of query gene'; grep "+GENE_ID+" "+GENE_BED+" > "+QUERY_BED+";"
    cmd_fastaFromBed="fastaFromBed -name -fi "+param['genome_fa']+" -bed "+QUERY_BED+" -fo "+QUERY_FA
    cmd=cmd+" "+cmd_fastaFromBed
    return(cmd)

def cmd_paired_end_for_query(param,query):
    inter_outdir = param['file_prefix']+query+'/'
    result_outdir = param['result_folder_prefix']+query+'/'
    log_outdir = param['LOG_F']+query+'/'
    
    #run_log?????? preprocess (it is not at this level), pair, single, summary each have one run log? And a log for each of these big step?
    log_outdir_paired_end = log_outdir+"log_of_paired_end_run/"
    log_outdir_single_end = log_outdir+"log_of_singleton_unmapped_run/"
    log_outdir_summary = log_outdir+"log_of_summary_run/"
    
    GENE_BED=param['whole_gene_list']
    QUERY_BED=result_outdir+"query_gene.bed"
    QUERY_FA=result_outdir+"query_gene.fa"    
    
    QF_pair_end_process_cmd="python "+param["QF_path"]+"QF_pair_end_process.py"+" -i "+inter_outdir+" -B "+param["bam_dir"]+" -o "+result_outdir+" -w "+param["whole_gene_list"]+" -g "+log_outdir_paired_end+" -t "+param["tophat_genome_ref"]+" -T "+param["genome_fa"]+" -F "+param["QF_path"]+" -l "+str(param["read_len"])+" -r "+str(param["resume_stat"])+" -a "+str(param["Align_percent"])+" -Q "+QUERY_BED+" -q "+param["step_size_query"]
    #print QF_pair_end_process_cmd
    return(QF_pair_end_process_cmd)

def cmd_single_end_for_query(param,query):
    inter_outdir = param['file_prefix']+query+'/'
    result_outdir = param['result_folder_prefix']+query+'/'
    log_outdir = param['LOG_F']+query+'/'
    
    #run_log?????? preprocess (it is not at this level), pair, single, summary each have one run log? And a log for each of these big step?
    log_outdir_paired_end = log_outdir+"log_of_paired_end_run/"
    log_outdir_single_end = log_outdir+"log_of_singleton_unmapped_run/"
    log_outdir_summary = log_outdir+"log_of_summary_run/"
    
    GENE_BED=param['whole_gene_list']
    QUERY_BED=result_outdir+"query_gene.bed"
    QUERY_FA=result_outdir+"query_gene.fa"
    
    QF_single_end_process_cmd="python "+param["QF_path"]+"QF_single_end_process.py"+" -i "+inter_outdir+" -B "+param["bam_dir"]+" -o "+result_outdir+" -w "+param["whole_gene_list"]+" -g "+log_outdir_single_end+" -t "+param["tophat_genome_ref"]+" -T "+param["genome_fa"]+" -F "+param["QF_path"]+" -l "+str(param["read_len"])+" -r "+str(param["resume_stat"])+" -a "+str(param["Align_percent"])+" -Q "+QUERY_BED+" -U "+param["step_size_query"]+" -O "+param["step_size_other"]+" -q "+query
    
    #print QF_single_end_process_cmd
    return(QF_single_end_process_cmd)

def summary_as_whole(param):
    # Extract real values from param dict
    num_parallels = int(param['num_processors'])

    queries = param['stub']
    timeout = int(param['timeout'])
        
    processes_sum = {}
    results_sum = {}
    cur_num_procs = 0
    # Dispatch each query
    for query in queries:
        #setup parameters:
        #run create output directories (mainly will go into two divisions: intermediate and result.)
        BAM_dir = param['bam_dir']
        inter_dir = param['file_prefix']
        result_dir = param['result_folder_prefix']
        
        #BAM_outdir = param['bam_dir']+param['stub'][param['file_index']]+'/'
        inter_outdir = param['file_prefix']+query+'/'
        result_outdir = param['result_folder_prefix']+query+'/'
        log_outdir = param['LOG_F']+query+'/'
        
        #run_log?????? preprocess (it is not at this level), pair, single, summary each have one run log? And a log for each of these big step?
        log_outdir_paired_end = log_outdir+"log_of_paired_end_run/"
        log_outdir_single_end = log_outdir+"log_of_singleton_unmapped_run/"
        log_outdir_summary = log_outdir+"log_of_summary_run/"
        
        #if not os.path.exists(BAM_outdir):
        #    os.makedirs(BAM_outdir)
        #check path and create folders
        if not os.path.exists(inter_outdir):
            os.makedirs(inter_outdir)
        
        if not os.path.exists(result_outdir):
            os.makedirs(result_outdir)
            
        if not os.path.exists(log_outdir):
            os.makedirs(log_outdir)
        
        if not os.path.exists(log_outdir_paired_end):
            os.makedirs(log_outdir_paired_end)
        
        if not os.path.exists(log_outdir_single_end):
            os.makedirs(log_outdir_single_end)
           
        if not os.path.exists(log_outdir_summary):
            os.makedirs(log_outdir_summary)
        
        # If there are not enough free processes
        while int(cur_num_procs) + 1 > int(num_parallels):
            # Wait for enough processes to terminate
            print "Hold at the waiting line"
            cur_num_procs -= wait_for(processes_sum, results_sum, timeout)
        # Launch jobs for this query
        proc_summary = async_cmd(cmd_summary_end_for_query(param,query))
        # Store jobs about this query in a dictionary
        jobs = dict(summary = proc_summary)
        processes_sum[query] = jobs
        # Increment number of processes active
        cur_num_procs += 1
    # Wait for all extant jobs to finish
    while(cur_num_procs != 0):
        cur_num_procs -= wait_for(processes_sum, results_sum, timeout)
    #print(processes_sum)
    #print(results_sum)
    #for queries in results_sum.keys():
    #    print queries
    #    for jobs in results_sum[queries].keys():
    #        print jobs
    #        for logs in results_sum[queries][jobs]:
    #            if type(logs) is types.IntType:
    #                print logs
    #            else:
    #                print logs[0]

def cmd_summary_end_for_query(param,query):
    inter_outdir = param['file_prefix']+query+'/'
    result_outdir = param['result_folder_prefix']+query+'/'
    log_outdir = param['LOG_F']+query+'/'
    
    #run_log?????? preprocess (it is not at this level), pair, single, summary each have one run log? And a log for each of these big step?
    log_outdir_paired_end = log_outdir+"log_of_paired_end_run/"
    log_outdir_single_end = log_outdir+"log_of_singleton_unmapped_run/"
    log_outdir_summary = log_outdir+"log_of_summary_run/"
    
    GENE_BED=param['whole_gene_list']
    QUERY_BED=result_outdir+"query_gene.bed"
    QUERY_FA=result_outdir+"query_gene.fa"
    
    QF_summary_process_cmd="python "+param["QF_path"]+"QF_summary_process.py"+" -i "+inter_outdir+" -B "+param["bam_dir"]+" -o "+result_outdir+" -w "+param["whole_gene_list"]+" -g "+log_outdir_summary+" -t "+param["tophat_genome_ref"]+" -T "+param["genome_fa"]+" -F "+param["QF_path"]+" -l "+str(param["read_len"])+" -r "+str(param["resume_stat"])+" -a "+str(param["Align_percent"])+" -Q "+QUERY_BED+" -s "+str(param["read_std"])+" -L "+str(param["split_n"])+" -N "+str(param["span_n"])+" -u "+str(param["sum_n"])+" -f "+str(param["span_only_filter"])+" -U "+param["step_size_query"]+" -O "+param["step_size_other"]
    return(QF_summary_process_cmd)
        
# Queue up a job, and immediately stop and wait for its results
def sync_cmd(cmd):    
    proc = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    retcode = proc.returncode
    return (stdout, stderr, retcode)

# Queue up a job, and get an object describing it which we can later ask for results from
def async_cmd(cmd):
    proc = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return proc

#need to build a log integrate function for both sync_cmd and async_cmd
#def log_integrate():

# Busy Wait for a job to finish
def wait_for(processes, results, timeout = 0):
    print("---------- Begin Busy Waiting ----------")
    waiting = True
    timer = 0
    finished = []
    total_done = 0
    # Busy Waiting Loop
    while(waiting):
        #print("Polling all processes")
        #print processes
        # Poll all processes
        for query, procs in processes.items():
            is_done = all([(proc.poll() is not None) for proc in procs.values()])
            # If 
            if is_done:
                if not query in results.keys():
                    results[query]={}
                h_overall_log=open(param['LOG_F']+query+'/'+"whole_procedure.log","a")
                for job_type, job in procs.items():
                    results[query][job_type] = (job.communicate(), job.returncode)
                    stdout_now=results[query][job_type][0][0]
                    stderr_now=results[query][job_type][0][1]
                    returncode_now=results[query][job_type][1]
                    h_overall_log.write(job_type+"\n"+stdout_now+"\n"+stderr_now+"\n")
                    if returncode_now!=0:
                        QF_all_modules.writeLog(job_type+' process in :'+query+' failed.\n',param)
                        print(job_type+' process in :'+query+' failed in each query module.\n Please see log files for details.')
                    else:
                        QF_all_modules.writeLog('finished processing '+job_type+' in :'+query+'.\n',param)
                h_overall_log.close()
                QF_all_modules.writeLog(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+" \n\n",param)
                #results[query] = {[job_type:(job.communicate(), job.returncode) for job_type, job in procs]}
                total_done += len(results[query])
                print("Total Done: %d" % total_done)
                finished.append(query)
                waiting = False
        timer += 10
        if timeout > 0 and timer > timeout:
            raise Exception("Waited too long for jobs, %r, %r" % (processes, results))
        #sleep can control how often to check in every sec(s)
        sleep(10)
    for fin in finished:
        #to delete the finished jobs in processes list
        processes.pop(fin)
    print("----------End Busy Waiting----------")
    return total_done

if __name__ == "__main__":
## Import modules
    #import matplotlib
    #matplotlib.use('Agg')
    import os, sys, shutil, time, shlex, subprocess, random, getopt, types
    #import json
    import QF_all_modules
    from time import sleep
    #sys.path.insert(0, 'path')
    #import module
    
    list_args = sys.argv[1:]
    
    if '-h' in list_args:
        print __doc__
        sys.exit(1)
    elif len(list_args) <2 or '-p' not in list_args:
        print __doc__
        sys.exit(1)
    else:
        param, parameter_file, updates = QF_all_modules.update_parameters(list_args)
       
    
    #if we resume we specify the old parameter file
    QF_all_modules.initialize_standard(param)
    QF_all_modules.write_updated_file(updates, param, parameter_file)

    #QF_all_modules.initialize_qsub(param)#looks like this is not needed any more since I am not using nested qsub to submit jobs but using waiting for loop function.
    QF_all_modules.readQueryIDFilenames(param)
    QF_all_modules.initialize_logfiles(param)
    
    #writeLog('Initializing all module parameters ... \n',param)
    #this step is to check all the parameter, but I do not have that many parameters to check, should be easier.
    #initialize_all(param) 

    #writeLog('Initializing successful!\n',param)    
    #writeLog('####################################################\n',param)    
    
    if param["clean_run"]=="FALSE":
	param["resume_stat"]=1
    else:
	param["resume_stat"]=0
    
    #set up folders path
    if not os.path.exists(param['working_dir']):
        os.makedirs(param['working_dir'])
    param['bam_dir']=param['working_dir']+'bams/'
    if not os.path.exists(param['bam_dir']):
        os.makedirs(param['bam_dir'])
    param['file_prefix']=param['working_dir']+'intermedias/'
    if not os.path.exists(param['file_prefix']):
        os.makedirs(param['file_prefix'])
    param['result_folder_prefix']=param['working_dir']+'results/'
    if not os.path.exists(param['result_folder_prefix']):
        os.makedirs(param['result_folder_prefix'])
    param['LOG_F']=param['working_dir']+'logs/'
    if not os.path.exists(param['LOG_F']):
        os.makedirs(param['LOG_F'])
    param["QF_path"]=param["scripts_dir"]
    
## run preprocess
    preprocess_log=param['result_folder_prefix']+'Preprocess.log'
    QF_all_modules.writeLog('Running all modules: \n\n',param)
    step_name="preprocessing in QF_wrapper.py"
    QF_all_modules.writeLog('preprocessing in QF_wrapper.py \n',param)
    next_step_name="finished preprocessing in QF_wrapper.py"
    #preprocess_cmd="python "+param["QF_path"]+"QF_preprocess.py "+param["accepted_bam"]+" "+param["unmapped_bam"]+" "+param['bam_dir']+" "+param["Perl_path"]+" "+param["QF_path"]+" "+param["whole_gene_list"]
    preprocess_cmd="python "+param["QF_path"]+"QF_preprocess.py "+param["accepted_bam"]+" "+param["unmapped_bam"]+" "+param['bam_dir']+" "+" "+param["QF_path"]+" "+param["whole_gene_list"]
    #if the bam files are splited, do not do it again.
    preprocess_status=QF_all_modules.resume_func(preprocess_cmd, 1, step_name, next_step_name, param['log_handle'])
    #check whether this step failed, if yes, stop the program.
    if preprocess_status[3]==0:
        QF_all_modules.writeLog('preprocessing in QF_wrapper.py done.\n',param)
        QF_all_modules.writeLog(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+" \n",param)
        h_preprolog_OUT=open(preprocess_log,"w")
        h_preprolog_OUT.write(preprocess_status[1])
        h_preprolog_OUT.close()
    else:
        h_preprolog_OUT=open(preprocess_log,"w")
        h_preprolog_OUT.write(preprocess_status[1]+preprocess_status[2])
        h_preprolog_OUT.close()
        QF_all_modules.writeLog('preprocessing failed, please see the results/Preprocess.log',param)
        sys.exit('Preprocess Failed, please see the results/Preprocess.log')

## run for each single query
    #run queryfuse for each query without summary
    each_query(param)
    
    #run summary for each query in multi node mode.
    if param['num_processors']>1:
        summary_as_whole(param)
    
    #need to add report for all querys
    
    #report all the queries together
    #pipeline.report_all(param) need to think more about how to report.
    #helper.report_finish(param)
    
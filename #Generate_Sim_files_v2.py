#Generate_Sim_files.py
import os, sys
import logging

def writeNewPreSolver(old_svpre_name,new_svpre_name,KW,):
	old = open(old_svpre_name, 'r')
	new = open(new_svpre_name,'w')
	
	line = old.readline()
	while line:
		if(line.startswith('set_surface_id_vtp mesh-complete/mesh-surfaces/### 199')):
			line = 'set_surface_id_vtp mesh-complete/mesh-surfaces/'+KW+' 199'
			line += '\n'
		if(line.startswith('pressure_vtp mesh-complete/mesh-surfaces/### 0')):
			line = 'pressure_vtp mesh-complete/mesh-surfaces/'+KW+' 0'
			line += '\n'
		new.write(line)
		line = old.readline()
	old.close()
	new.close()

def writeNewBatch(old_name,new_name,KW,):
	old = open(old_name, 'r')
	new = open(new_name,'w')
	
	line = old.readline()
	while line:
		if(line.find('job-name')!=-1):
			line = '#SBATCH --job-name='+KW
			line += '\n'
		new.write(line)
		line = old.readline()
	old.close()
	new.close()

cwd = os.getcwd()
parentdir = os.path.join(cwd,'Collateral_Sims')
Sim_results = os.path.join(parentdir,'Collateral_results')
if(not os.path.exists(parentdir)):
	logging.info('Making %s ',parentdir)
	os.mkdir(parentdir)
if(not os.path.exists(Sim_results)):
	logging.info('Making %s ',Sim_results)
	os.mkdir(Sim_results)
logfile = os.path.join(parentdir,'Simulations_logfile.log')
logging.basicConfig(filename=logfile,level=logging.DEBUG)
flow_file = os.path.join(cwd,'aorta.flow')
svpre_file = os.path.join(cwd,'Sim.svpre')
solver_file = os.path.join(cwd,'solver.inp')
cort_file = os.path.join(cwd,'cort.dat')
rcr_file = os.path.join(cwd,'rcrt.dat')
presolver = '~/SimVascular3.0.170301/Bin/svpre'
presolver_comet = '~/svsolver/BuildWithMake/Bin/svpre.exe'
postsolver_comet = '~/svsolver/BuildWithMake/Bin/svpost.exe'
batch_file = os.path.join(cwd,'flowsolver')
batch_file_comet = '/oasis/scratch/comet/suhaasa/temp_project/flowsolver'
procs_folder = '96-procs_case'
scale_script = 'scale_script.sh'
presolve = 1
solve = 0
postsolve = 0

j_strings = ['lca','lca75','lca80','lca85','lca90','lca95','lca97','lca99']
j_strings_svpre = ['LCA_main','LCA_main_75','LCA_main_80','LCA_main_85','LCA_main_90','LCA_main_95','LCA_main_97','LCA_main_99']
i_strings = ['0col','3cola_40um','6col_20um','3colb_40um','12col_20um','6col_28um','3col_40um','16col_20um','1col5_40um','1col9_40um','1col10_40um']
Mesh_set = set()
PreSim_set = set()
Sim_set = set()
PostSim_set = set()
k=0
tag = ['','a','b']
for i in range(0,len(i_strings)):
	for j in range(0,len(j_strings)):
		Sim_code = 'Sim'+str(i)+str(j)
		Sim_name = 'Sim'+str(i)+str(j)+'_'+j_strings[j]+'_'+i_strings[i]
		meshdir = os.path.join('..','..','Meshes','Latest_meshes','full_model'+'_'+j_strings[j]+'_'+i_strings[i]+'-mesh-complete')
		Sim_path = os.path.normcase(os.path.join(parentdir,Sim_name))
		if(not os.path.exists(Sim_path)):
			logging.info('Making %s ',Sim_path)
			os.mkdir(Sim_path)
		else: 
			logging.info('%s exists.',Sim_path)
		logging.info('Copying %s into %s',flow_file,Sim_path)
		os.system('cp '+flow_file+' '+Sim_path)
		logging.info('Copying %s into %s',solver_file,Sim_path)
		os.system('cp '+solver_file+' '+Sim_path)
		os.system('cp '+rcr_file+' '+Sim_path)
		os.system('cp '+cort_file+' '+Sim_path)
		if(not os.path.exists(os.path.join(Sim_path,Sim_code+'.svpre'))):
			logging.info('Making %s ',os.path.join(Sim_path,Sim_code+'.svpre'))
			writeNewPreSolver(svpre_file,os.path.join(Sim_path,Sim_code+'.svpre'),j_strings_svpre[j]+'.vtp')
		else:
			logging.info('%s exists.',os.path.join(Sim_path,Sim_code+'.svpre'))
		
		if(os.path.exists(meshdir) and not os.path.exists(os.path.normcase(os.path.join(Sim_path,'mesh-complete')))):
			logging.info('Copying %s into %s and renaming to %s',meshdir,Sim_path,'mesh-complete')
			os.system('cp -r '+meshdir+' '+Sim_path)
			os.system('mv '+os.path.normcase(os.path.join(Sim_path,os.path.basename(meshdir)))+' '+os.path.normcase(os.path.join(Sim_path,'mesh-complete')))
			os.system('mv '+os.path.normcase(os.path.join(Sim_path,'mesh-complete','full*.vtu'))+' '+os.path.normcase(os.path.join(Sim_path,'mesh-complete','mesh-complete.mesh.vtu')))
			os.system('mv '+os.path.normcase(os.path.join(Sim_path,'mesh-complete','full*.vtp'))+' '+os.path.normcase(os.path.join(Sim_path,'mesh-complete','mesh-complete.exterior.vtp')))
			os.system('cp '+os.path.normcase(os.path.join(cwd,scale_script))+' '+Sim_path)
			os.chdir(Sim_path)
			os.system('sh '+scale_script+' 1.25')	
		if(os.path.exists(os.path.normcase(os.path.join(Sim_path,'mesh-complete')))):
			logging.info('%s has a mesh-complete dir.',Sim_name)
			Mesh_set.add(Sim_name)
		if(presolve and Sim_name in Mesh_set and not os.path.exists(os.path.join(Sim_path,'restart.0.1'))):
			logging.info('Changing dir to %s',Sim_path)
			os.chdir(Sim_path)
			print('Running presolver on '+Sim_code)
			logging.info('Running presolver on %s',Sim_code)
			os.system(presolver+' '+Sim_code+'.svpre')
		if(solve):
			os.chdir(Sim_path)
			#if(os.path.exists(os.path.normcase(os.path.join(Sim_path,procs_folder)))):
			#	os.system('mv '+procs_folder+'_old'+' '+procs_folder)
			writeNewBatch(batch_file,'./flowsolver_'+Sim_code,Sim_code)
			os.system('sbatch flowsolver_'+Sim_code)
		if(os.path.exists(os.path.join(Sim_path,'restart.0.1'))):
			logging.info('Presolve successfully completed for %s',Sim_code)
			PreSim_set.add(Sim_name)

		if(postsolve):
			logging.info('Running postsolver on %s',Sim_code)
			os.system(postsolver_comet+' -indir '+os.path.join(Sim_path,procs_folder)+ ' -outdir '+os.path.join(Sim_results)+' -start 3000 -stop 4000 -incr 50 -vtu all_results_'+Sim_name+' -vtp all_results_'+Sim_name)
		if(os.path.exists(os.path.join(Sim_results,'all_results_'+Sim_name+'.vtu'))):
			logging.info('Results successfully stored for %s in %s',Sim_name,Sim_results)

		logging.info('Changing dir to %s',cwd)
		os.chdir(cwd)
print(Sim_set)
logging.info('Mesh ready: ' + ','.join(sorted(Mesh_set)))
logging.info('Simulation ready: ' + ','.join(sorted(PreSim_set)))
logging.info('Simulation done: ' + ','.join(sorted(Sim_set)))
logging.info('PostSimulation done: ' + ','.join(sorted(PostSim_set)))
print('done.')

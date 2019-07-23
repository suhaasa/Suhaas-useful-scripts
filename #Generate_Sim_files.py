#Generate_Sim_files.py
import os, sys
import logging

def writeNewPreSolver(old_svpre_name,new_svpre_name,KW,):
	old = open(old_svpre_name, 'r')
	new = open(new_svpre_name,'w')
	
	line = old.readline()
	while line:
		if(line.startswith('set_surface_id_vtp mesh-complete/mesh-surfaces/### 196')):
			line = 'set_surface_id_vtp mesh-complete/mesh-surfaces/'+KW+' 196'
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
logfile = os.path.join(parentdir,'Simulations_logfile.log')
logging.basicConfig(filename=logfile,level=logging.DEBUG)
flow_file = os.path.join(cwd,'aorta.flow')
svpre_file = os.path.join(cwd,'Sim.svpre')
solver_file = os.path.join(cwd,'solver.inp')
presolver = '~/SimVascular3.0.170301/Bin/svpre'
presolver_comet = '~/svsolver/BuildWithMake/Bin/svpre.exe'
postsolver_comet = '~/svsolver/BuildWithMake/Bin/svpost.exe'
batch_file = '/oasis/scratch/comet/suhaasa/temp_project/flowsolver'
procs_folder = '96-procs_case'
presolve = 0
postsolve = 0
if(not os.path.exists(parentdir)):
	logging.info('Making %s ',parentdir)
	os.mkdir(parentdir)
if(not os.path.exists(Sim_results)):
	logging.info('Making %s ',Sim_results)
	os.mkdir(Sim_results)
j_strings = ['lca','lca75','lca80','lca85','lca90','lca95','lca97','lca99']
j_strings_svpre = ['LCA_main','LCA_main_75','LCA_main_80','LCA_main_85','LCA_main_90','LCA_main_95','LCA_main_97','LCA_main_99']
i_strings = ['0col','3col_20um','6col_20um','9col_20um','12col_20um','6col_28um','3col_40um']
Mesh_set = set()
PreSim_set = set()
Sim_set = set()
PostSim_set = set()
for i in range(0,len(i_strings)):
	for j in range(0,len(j_strings)):
		Sim_name = 'Sim'+str(i)+str(j)+'_'+j_strings[j]+'_'+i_strings[i]
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
		if(not os.path.exists(os.path.join(Sim_path,'Sim'+str(i)+str(j)+'.svpre'))):
			logging.info('Making %s ',os.path.join(Sim_path,'Sim'+str(i)+str(j)+'.svpre'))
			writeNewPreSolver(svpre_file,os.path.join(Sim_path,'Sim'+str(i)+str(j)+'.svpre'),j_strings_svpre[j]+'.vtp')
		else:
			logging.info('%s exists.',os.path.join(Sim_path,'Sim'+str(i)+str(j)+'.svpre'))
		meshdir = os.path.join('..','Meshes','full_model'+'_'+j_strings[j]+'_'+i_strings[i]+'-mesh-complete')
		if(os.path.exists(meshdir) and not os.path.exists(os.path.normcase(os.path.join(Sim_path,'mesh-complete')))):
			logging.info('Copying %s into %s and renaming to %s',meshdir,Sim_path,'mesh-complete')
			os.system('cp -r '+meshdir+' '+Sim_path)
			os.system('mv '+os.path.normcase(os.path.join(Sim_path,os.path.basename(meshdir)))+' '+os.path.normcase(os.path.join(Sim_path,'mesh-complete')))
		if(os.path.exists(os.path.normcase(os.path.join(Sim_path,'mesh-complete')))):
			logging.info('%s has a mesh-complete dir.',Sim_name)			
			Mesh_set.add(Sim_name)
		if(presolve and not os.path.exists(os.path.join(Sim_path,'restart.0.1')) and not os.path.exists(meshdir) and os.path.exists(os.path.normcase(os.path.join(Sim_path,'mesh-complete')))):
			logging.info('Changing dir to %s',Sim_path)
			os.chdir(Sim_path)
			print('Running presolver on '+'Sim'+str(i)+str(j))
			logging.info('Running presolver on %s','Sim'+str(i)+str(j))
			os.system(presolver_comet+' '+'Sim'+str(i)+str(j)+'.svpre')
			if(os.path.exists(os.path.join(Sim_path,'restart.0.1'))):
				writeNewBatch(batch_file,'./flowsolver_'+'Sim'+str(i)+str(j),'Sim'+str(i)+str(j))
				os.system('sbatch flowsolver_'+'Sim'+str(i)+str(j))
		if(os.path.exists(os.path.join(Sim_path,'restart.0.1'))):
			logging.info('Presolve successfully completed for %s','Sim'+str(i)+str(j))
			PreSim_set.add(Sim_name)

		if(postsolve and os.path.exists(os.path.join(Sim_path,procs_folder))):
			logging.info('Running postsolver on %s','Sim'+str(i)+str(j))
			os.system(postsolver_comet+' -indir '+os.path.join(Sim_path,procs_folder)+ ' -outdir '+os.path.join(Sim_results)+' -start 100 -stop 100 -incr 0 -vtu all_results_'+Sim_name)
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
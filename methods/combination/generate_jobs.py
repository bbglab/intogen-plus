import os
cancer_types = ["ACC","BLCA","BRCA","CESC","CHOL","COAD","COADREAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]
print ("#source activate Schulze")
dir_output = "/workspace/projects/intogen/intogen4/scripts/data/results/optimization/ranking/"
seeds = "12"
niter = "1000"
epsilon = "0.1"
input_rankings = "/workspace/projects/intogen/intogen4/scripts/data/dict_parsed_methods_ranking.pickle"
discarded_methods = "/workspace/projects/intogen/intogen4/runs/intogen4_20170614/plots/discarted_analyses.txt"
name = "RANKING_WORSTCASE"
dir_input = "/workspace/projects/intogen/intogen4/scripts/data/results/optimization/weights/"
def generate_jobs_pan_cancer():

	for cancer in cancer_types:
		if  not os.path.exists(dir_output+"/"+cancer):
			os.mkdir(dir_output+"/"+cancer)
		print ("python /workspace/projects/intogen/intogen4/scripts/Schulze/optimizer.py --seeds "+seeds+" --niter "+niter+" --epsilon "+epsilon+" --foutput "+dir_output+cancer+"/weights --input_rankings "+input_rankings+" --cancer "+cancer + " --discarded_methods "+discarded_methods+" --t_combination "+"RANKING")

def generate_jobs_cv(N=20,percentage=0.8):

	for cancer in cancer_types:
		if not os.path.exists(dir_output+cancer):
			os.mkdir(dir_output+cancer)
		for i in range(N):
			print ("python /workspace/projects/intogen/intogen4/scripts/Schulze/optimizer.py --seeds "+seeds+" --niter "+niter+" --epsilon "+epsilon+" --foutput "+dir_output+"/"+cancer+"/"+str(i)+" --input_rankings "+input_rankings+" --cancer "+cancer + " --discarded_methods "+discarded_methods+" --t_combination "+name +"  --percentage_cgc "+"0.80")



def generate_jobs_schulze(N=98,name="ranking"):
	#for cancer in cancer_types:

	print ("python /workspace/projects/intogen/intogen4/scripts/Schulze/schulze.py --input_data "+input_rankings+" --directory_optimize_weights "+dir_input+" --directory_output  "+dir_output+" --dict_output "+dir_output+"dict_ranking_optimized.pickle"+" --name Combination_Ranking_Optimized --type_run optimization --type_input ranking")


def generate_random_voting():

	print ("python /workspace/projects/intogen/intogen4/scripts/Schulze/schulze.py --input_data "+input_rankings+" --directory_optimize_weights "+dir_input+" --directory_output  "+dir_output+" --dict_output "+dir_output+"ranking_default_weights.pickle"+" --name Combination_Ranking --type_run default --type_input ranking")
generate_jobs_schulze()
#generate_jobs_cv()
#generate_random_voting()
#generate_jobs_pan_cancer()

import pandas as pd
import numpy as np
import pickle
#from Parser import Parser
import Utils
class Evaluation_Simulated:


    @staticmethod
    def read_simulated(path_simulated):
        '''
        Read the simulated drivers and returns a DataFrame with the information of drivers per cancer type
        :param path_simulated: path of the simulated data
        :return: the DataFrame with the information
        '''
        df_simulated = pd.read_csv(path_simulated,sep="\t")
        df_simulated.rename(columns={"Unnamed: 0":"Cancer_Type"},inplace=True)
        df_simulated["N_Drivers"] = df_simulated.apply(lambda row: int(row["Simulated drivers"].split(" ")[1]),axis=1)
        df_simulated["Drivers"] = df_simulated.apply(lambda row: row["Simulated drivers"].split(" ")[2][1:-1],axis=1)
        return df_simulated


    def __init__(self,path_simulated="/workspace/projects/intogen/intogen4/report_simulations.txt"):

        self.df_simulated = Evaluation_Simulated.read_simulated(path_simulated)
        self.type_evaluation = "Evaluation_Simulated"


    def get_positives_negatives(self,list_ranking,cancer,list_total):
        '''
        Given a list of all predictions returns the TP,FP,FN,TN
        :param list_ranking: list of ranked genes (P)
        :param cancer: Cancer Type
        :param list_total: list of all predictions by the method (include P and N)
        :return: the TP,FP,FN,TN
        '''
        if self.df_simulated[self.df_simulated["Cancer_Type"]==cancer]["N_Drivers"].values[0] == 0:
            tn = len(list_total) - (len(list_ranking))
            return 0,len(list_ranking),0,tn # tp,fp,fn,tn
        total_tp = self.df_simulated[self.df_simulated["Cancer_Type"]==cancer]["Drivers"].values[0].split(",")
        tpl = [gene for gene in list_ranking if gene in total_tp]
        tp = len(tpl)
        fp = len(list_ranking) - tp
        fn = self.df_simulated[self.df_simulated["Cancer_Type"]==cancer]["N_Drivers"].values[0] - tp
        tn = len(list_total) - (tp+fp+fn)
        return float(tp),float(fp),float(fn),float(tn)


    def get_precision_recall(self,tp,fp,fn,tn):
        '''
        Given the TP,FP,FN,TN returns the Precision,Recall,F-Measure and MCC
        :param tp: True positives i.e. genes hit by the method that are drivers
        :param fp: False positives i.e. genes hit by the method that had not been simulated as drivers
        :param fn: False negatives i.e. genes that are drivers and are not found by the method
        :param tn: True negatives i.e. genes that are not simulated as drivers and are not found in the predictions of the method.
        :return: Precision, Recall, F-Measure, MCC
        '''

        if tp + fp ==0 or tp+fn == 0:
            precision = 0.0
            recall = 0.0
            f1 = 0.0
            MCC = 0.0
        else:
            precision = float(tp / (tp+fp))
            recall = float(tp / (tp+fn))
            f1 = 2*tp / (2*tp + fn + fp)
            MCC = ((tp*tn) - (fp*fn)) / np.sqrt((tp+fn)*(tp+fp)*(tn+fp)*(tn+fn))
        return precision,recall,f1,MCC

    l_stats = []
    def get_list_ranking(self,d_results):
        '''
        Given a dictionary with the ranking of the genes returns a list with keeping the ranking of the genes.
        :param d_results: dictionary of format {gene:ranking,gene2:ranking2,gene3:ranking3....}
        :return: a list with the sorted ranking of genes
        '''

        return [x[0] for x in sorted(d_results.items(),key=lambda x: (x[1],x[0]))]

    def calculate_stats_methods_global(self,dict_results,dict_total,methods= ["hotmapssignature","oncodrivefml","mutsigcv","oncodriveomega","oncodriveclust"],cancer_types = ["ACC","BLCA","BRCA","CESC","CHOL","COAD","COADREAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]):
        '''
        Calculate the stats globally for the input dict_results
        :param dict_results: dict of results with the format {cancer:{method:{"gene":ranking,"gene2":ranking...}}
        :param dict_total: dict with all the predictions with the format {cancer:{method:{"gene":ranking,"gene2":ranking...}}
        :param methods: list of query methods
        :param cancer_types: list of query cancer types
        :return: return a DataFrame with the stats of the the evaluation
        '''
        l_stats = []
        for cancer in cancer_types:

            for method in methods:
                if method in dict_results[cancer]:

                    list_query = self.get_list_ranking(dict_results[cancer][method])
                    list_total = self.get_list_ranking(dict_total[cancer][method])
                    if len(list_total) < len(list_query):
                        list_total = list_query
                    tp,fp,fn,tn = self.get_positives_negatives(list_query,cancer,list_total)

                    precision,recall,f1,MCC = self.get_precision_recall(tp,fp,fn,tn)
                    l_stats.append([cancer,method,precision,recall,f1,MCC,tp,fp,fn,tn])


        df_stats = pd.DataFrame(l_stats,columns=["Cancer_Type","Method","Precision","Recall","F1","MCC","TP","FP","FN","TN"])
        return df_stats


    def calculate_stats_methods_step(self,dict_results,dict_total,begin=5,end=60,step=5,methods= ["hotmapssignature_t","oncodrivefml_t","mutsigcv_t","oncodriveomega_t","oncodriveclust_t","hotmapssignature_r","oncodrivefml_r","mutsigcv_r","oncodriveomega_r","oncodriveclust_r","Combination_Threshold","Combination_Ranking"],cancer_types = ["ACC","BLCA","BRCA","CESC","CHOL","COAD","COADREAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]):
        '''
        Calculate the stats for different step sizes of the input dictionaries.
        :param dict_results: dict of results with the format {cancer:{method:{"gene":ranking,"gene2":ranking...}}
        :param dict_total: dict with all the predictions with the format {cancer:{method:{"gene":ranking,"gene2":ranking...}}
        :param methods: list of query methods
        :param cancer_types: list of query cancer types
        :return: return a DataFrame with the stats of the the evaluation
        '''
        l_stats = []
        for i in range(begin,end,step):

            for cancer in cancer_types:

                for method in methods:
                    if method in dict_results[cancer]:

                        l= self.get_list_ranking(dict_results[cancer][method])
                        list_query = list(l[0:i])
                        if not method in dict_total[cancer]:
                            list_total = list(l)
                        else:
                            list_total = self.get_list_ranking(dict_total[cancer][method])

                        tp,fp,fn,tn = self.get_positives_negatives(list_query,cancer,list_total)

                        precision,recall,f1,MCC = self.get_precision_recall(tp,fp,fn,tn)
                        l_stats.append([cancer,method,precision,recall,f1,MCC,tp,fp,fn,tn,i])


        df_stats = pd.DataFrame(l_stats,columns=["Cancer_Type","Method","Precision","Recall","F1","MCC","TP","FP","FN","TN","TOP"])
        return df_stats

if __name__ == "__main__":
    d_results_methodsr= pickle.load( open( "/workspace/projects/intogen/intogen4/scripts/data/dict_parsed_methods_ranking_simulated.pickle", "rb" ) )
    d_results_methodst= pickle.load( open( "/workspace/projects/intogen/intogen4/scripts/data/dict_parsed_methods_threshold_simulated.pickle", "rb" ) )
    d_results_cranking = pickle.load( open( "/workspace/projects/intogen/intogen4/scripts/data/ranking_combination_simulated.pickle", "rb" ) )
    d_results_cthreshold = pickle.load( open( "/workspace/projects/intogen/intogen4/scripts/data/threshold_combination_simulated.pickle", "rb" ) )
    d_results = Utils.join_dictionaries([d_results_methodsr,d_results_methodst,  d_results_cranking,d_results_cthreshold])
    d_total = pickle.load( open( "/workspace/projects/intogen/intogen4/scripts/data/dict_parsed_methods_all_simulated.pickle", "rb" ) )
    e = Evaluation_Simulated()
    print (d_results)
    df_stats = e.calculate_stats_methods_step(d_results,d_total)
    df_stats.to_csv("/workspace/projects/intogen/intogen4/scripts/data/stats_simulated.csv",sep="\t",index=False)



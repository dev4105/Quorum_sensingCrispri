import pandas as pd

class convert2df():
    
    def __init__(self):
        self.df = pd.DataFrame()
    
    def quorum(self, sol_array):
        temp_dataframe = pd.DataFrame(columns = ['inAcP', 'AcP', 'N', 'LuxR', 'AHL', 'LuxR_A', 'GFP']) 
        for i in range(len(temp_dataframe.columns)):
            temp_dataframe[temp_dataframe.columns[i]] = sol_array[:,i]
        self.df = temp_dataframe

    def qcrispri(self, sol_array):
        temp_dataframe = pd.DataFrame(columns = ['inAcP', 'AcP', 'N', 'LuxR', 'AHL', 'LuxR_A', 'GFPc', 'GFP', 'dCas9', 'pC', 'pD'])
        for i in range(len(temp_dataframe.columns)):
            temp_dataframe[temp_dataframe.columns[i]] = sol_array[:,i]
        self.df = temp_dataframe
    
    def qcrispriSimple(self, sol_array):
        temp_dataframe = pd.DataFrame(columns = ['AHL', 'dCas9', 'GFP', 'N'])
        for i in range(len(temp_dataframe.columns)):
            temp_dataframe[temp_dataframe.columns[i]] = sol_array[:,i]
        self.df = temp_dataframe
    
    
    def qcrispri_sponge(self, sol_array):
        temp_dataframe = pd.DataFrame(columns = ['inAcP', 'AcP', 'N', 'LuxR', 'AHL', 'LuxR_A', 'GFPc', 'GFP', 'dCas9', 'pC', 'pD', 'pS', 'pSB'])
        for i in range(len(temp_dataframe.columns)):
            temp_dataframe[temp_dataframe.columns[i]] = sol_array[:,i]
        self.df = temp_dataframe
    
    def returndf(self):
        return(self.df)
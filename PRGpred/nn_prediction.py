from keras.models import load_model, Sequential
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
import pkg_resources
from PRGpred import utils

def pred_Phase1(df, outfile):
    testData=df['Samples']
    Lables=df['Labels']
    SampleID=df['SeqID']
    Samples=testData.reshape(testData.shape[0], 1, testData.shape[1],1)
    file=pkg_resources.resource_filename('PRGpred','data/rgene_vs_non-rgene.h5')
    myModel=load_model(file,compile=False)
    prob=myModel.predict(Samples)
    pred=np.argmax(prob, axis=1)
    df = pd.DataFrame()
    df = df.assign(SampleID=pd.Series(SampleID))
    df = df.assign(prob_Rgenes=pd.Series(prob.T[0]))
    df = df.assign(prob_Non_Rgenes=pd.Series(prob.T[1]))
    df = df.assign(pred=pd.Series(pred))
    df['Prediction'] = np.where(df['pred'] == 0, 'Rgenes', 'Non-Rgenes')
    df.columns = ['SampleID', 'Rgenes', 'Non-Rgenes','pred', 'Prediction']
    dff = df[['SampleID', 'Prediction', 'Rgenes', 'Non-Rgenes']]
    s = round(dff.select_dtypes(include=[np.number]) * 100, 4)
    dff[s.columns] = s
    dff.to_csv(outfile,sep="\t", index=False)
    df2 = df.loc[df['Prediction'] =='Rgenes']
    dd= df2['SampleID'].values.tolist()

    return dd, dff


def pred_Phase2(df,outfile):
    """
    Function for Phase 1 prediction
    :param df:
    :return: Log file of Phase 1
    """
    testData=df['Samples']
    Lables=df['Labels']
    SampleID=df['SeqID']
    Samples=testData.reshape(testData.shape[0], 1, testData.shape[1],1)
    file=pkg_resources.resource_filename('PRGpred','data/rgene_classification.h5')
    myModel=load_model(file,compile=False)
    prob=myModel.predict(Samples)
    pred=np.argmax(prob, axis=1)
    df = pd.DataFrame()
    df = df.assign(SampleID=pd.Series(SampleID))
    df = df.assign(prob_Class1=pd.Series(prob.T[0]))
    df = df.assign(prob_Class2=pd.Series(prob.T[1]))
    df = df.assign(prob_Class3=pd.Series(prob.T[2]))
    df = df.assign(prob_Class4=pd.Series(prob.T[3]))
    df = df.assign(prob_Class5=pd.Series(prob.T[4]))
    df = df.assign(prob_Class6=pd.Series(prob.T[5]))
    df = df.assign(prob_Class7=pd.Series(prob.T[6]))
    df = df.assign(prob_Class8=pd.Series(prob.T[7]))
    df = df.assign(pred=pd.Series(pred))
    p3 = {0 : "CNL", 1: "KIN", 2: "LYK", 3: "LECRK", 
            4: "RLK", 5: "RLP",
            6:"TIR",7:"TNL"}
    pp= list(df['pred'])
    df['Prediction']=df['pred']
    df['Prediction'].replace(p3,inplace=True)
    
    enz = ['SampleID', 'CNL', 'KIN', 'LYK', 'LECRK','RLK', 'RLP', 'TIR', 'TNL', 'pred','Prediction']
    df.columns = enz
    s = round(df.select_dtypes(include=[np.number]) * 100,4)
    df[s.columns] = s
    df=df.drop('pred',1)
    df['pred']=pp
    dd=df[['SampleID','Prediction', 'CNL', 'KIN', 'LYK', 'LECRK','RLK', 'RLP', 'TIR', 'TNL']]
    dd.to_csv(outfile,sep="\t", index=False)
   
    
    return dd
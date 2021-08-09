# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 14:26:00 2021

@author: 19193
"""
import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import networkx as nx


data= pd.read_csv('Kinome_SIGNAL_(MinBknd).csv',header=None)
peptic_name=data.iloc[11:,1]
patient_id=np.unique(data.iloc[5,7:])
Array=np.unique(data.iloc[2,7:])
Cycle=np.unique(data.iloc[3,7:])
Exposure_time=np.unique(data.iloc[4,7:])
Drug=[ 'CT', 'Brigatinib', 'Sitravatinib','Neratinib']


sample_id=data.iloc[5,7:]+'_'+data.iloc[2,7:] +'_'+data.iloc[3,7:].astype('str')+'_'+data.iloc[4,7:].astype('str')

kinome_data=pd.DataFrame(data.iloc[11:,7:])
kinome_data.columns=sample_id
kinome_data.index=peptic_name

idx=np.where(data.iloc[2,:]=='A1')

A1= np.zeros((21, 194 ,68))

itr=0
for i in range(0+7,(68*4*21)+7,68*4):
    A1[itr,:,:]= data.iloc[11:,i:i+68]
    itr=itr+1
    
A11=np.median(A1,axis=0)

correlation=pd.DataFrame(A11.T).corr(method='pearson').abs()

correlation1=correlation[correlation>.99]

correlation1[correlation1.isnull()]=0

correlation1=np.array(correlation1)

correlation1=correlation1.flatten()

network= pd.DataFrame(np.zeros(( 194*194 ,3)))

network.iloc[:,2]=correlation1

peptic_name1=[]
for _ in range(194):
      peptic_name1.append(peptic_name)


peptic_name1=np.array(peptic_name1)

peptic_name1=peptic_name1.flatten()


network.iloc[:,0]=peptic_name1


peptic_name2=[]
for name in peptic_name:
    for _ in range(194):
      peptic_name2.append(name)

peptic_name2=np.array(peptic_name2)

network.iloc[:,1]=peptic_name2

idx=np.where(network.iloc[:,2]>0)[0]

network=network.iloc[idx,:]

idx=np.where(network.iloc[:,2]<1)[0]
network=network.iloc[idx,:]

network.columns =['Peptide_X', 'Peptide_Y', 'Sigma'] 
network = np.core.records.fromarrays(network.values.T,
                                          names='Peptide_X, Peptide_Y, Sigma',
                                          formats = 'S20, S20, f8')

G=nx.Graph()
G.add_weighted_edges_from(network)

pos= nx.kamada_kawai_layout(G)

G.nodes()

res=400
node_x = []
node_y = []
for node in G.nodes():
    x, y = pos[node]
    node_x.append(200*x+200)
    node_y.append(200*y+200)
nme=[]
for node in G.nodes():
    nme.append(node)


# sigma=.04

# sample_id=sample_id.to_list()

def cr_gt(data1,sigma,res,max_dat):  # Create for Gene terrain 
    
    x_ = np.linspace(-1, 1, res)
    y_ = np.linspace(-1, 1, res)
    X, Y = np.meshgrid(x_, y_)
    
    gaussian=np.zeros((res,res))
    
    pi=3.14
    itr=0
    for node in G.nodes():
        x, y = pos[node]
        amp1=np.exp(data1[itr]/max_dat)
        itr=itr+1
        gaussian1 =(1/np.sqrt(2*pi*sigma))*np.exp(-(((X-x)/sigma)**2+((Y-y)/sigma)**2)) 
        gaussian =amp1* gaussian1 + gaussian
        # gaussian = np.flipud(gaussian)
        gray=(gaussian-np.min(gaussian))/(np.max(gaussian)-np.min(gaussian))
    
    return 1-gray

st.title('KinomeTerrain')

st.header('Drug Sensitivity')

col1,col2,col3= st.beta_columns(3)
with col1:
    Patinet1 = st.selectbox(
                "Tumor", patient_id
            )

with col2:

    Cycle1 = st.selectbox(
                "Cycle", Cycle
            )
with col3:

    Exposure_time1 = st.selectbox(
                "Exposure_time", Exposure_time
            )

st.sidebar.header('Properties of Terrain')
sigma=st.sidebar.slider('Sigma',0.0,0.1,.04,.01)
btn= st.sidebar.checkbox('W/Peptide names')
res=st.sidebar.slider('Resulotion of image',100, 1000, 400,100)

cols = st.beta_columns(4) 

sample_id1=sample_id.to_list()
itr=0
max_dat=[]
for j in range(4):
            print(Patinet1+'_'+Array[j]+'_'+Cycle1+'_'+Exposure_time1)
            idx=np.where(sample_id==Patinet1+'_'+Array[j]+'_'+Cycle1+'_'+Exposure_time1)[0][0]
            print(idx)
            data1=kinome_data[sample_id1[idx]].astype('float')
            max_dat.append(np.max(data1))
            
            
print(max(max_dat))
max_dat=max(max_dat)


itr=0
for j in range(2):
        cols = st.beta_columns(2) 
        # re.columns=['MSE','SSIM']
        for i in range(2):
            print(Patinet1+'_'+Array[j]+'_'+Cycle1+'_'+Exposure_time1)
            idx=np.where(sample_id==Patinet1+'_'+Array[itr]+'_'+Cycle1+'_'+Exposure_time1)[0][0]
            data1=kinome_data[sample_id1[idx]].astype('float')
            gray = cr_gt(data1,sigma,res,max_dat)
            # [150:230,150:230]
            fig = px.imshow(gray,color_continuous_scale='spectral',width=400, height=400)
            if btn:
                fig.add_trace(go.Scatter(x=node_x, y=node_y,
                          text=nme,
                          textposition='top center',
                mode='text+markers',
                name='Point'
                          ))


            fig.update_layout(coloraxis_showscale=False)
            fig.update_layout( title=Array[itr]+ ' ( ' +Drug[itr]+' )')
            cols[i].plotly_chart(fig)
            itr=itr+1
              


# node_trace = go.Scatter(x=[], y=[], mode='markers+text',
#                      text=G.nodes(),
#                      textposition='top',
#                      )
# for node in G.nodes():
#     x, y = pos[node]
#     node_trace['x'].append(x)
#     node_trace['y'].append(y)
   


# st.plotly_chart(node_trace)

st.header('Comaprison between 2 Tumors in same condition')

col1,col2,col3= st.beta_columns(3)
with col1:
    Array1 = st.selectbox(
                "Array", Array,key=1
            )

with col2:

    Cycle1 = st.selectbox(
                "Cycle", Cycle,key=1
            )
with col3:

    Exposure_time1 = st.selectbox(
                "Exposure_time", Exposure_time,key=1
            )
    
col1,col2= st.beta_columns(2)
with col1:
    Patinet11 = st.selectbox(
                "Tumor 1", patient_id,key=1
            )
with col2:
    Patinet12 = st.selectbox(
                "Tumor 2", patient_id,key=2
            )

itr=0
max_dat=[]

print(Patinet1+'_'+Array[j]+'_'+Cycle1+'_'+Exposure_time1)
idx=np.where(sample_id==Patinet11+'_'+Array1+'_'+Cycle1+'_'+Exposure_time1)[0][0]

print(idx)
data1=kinome_data[sample_id1[idx]].astype('float')
max_dat.append(np.max(data1))
            

print(Patinet1+'_'+Array[j]+'_'+Cycle1+'_'+Exposure_time1)
idx=np.where(sample_id==Patinet12+'_'+Array1+'_'+Cycle1+'_'+Exposure_time1)[0][0]

print(idx)
data1=kinome_data[sample_id1[idx]].astype('float')
max_dat.append(np.max(data1))            

print(max(max_dat))
max_dat=max(max_dat)



itr=0



cols_p2 = st.beta_columns(2) 
idx=np.where(sample_id==Patinet11+'_'+Array1+'_'+Cycle1+'_'+Exposure_time1)[0][0]
data1=kinome_data[sample_id1[idx]].astype('float')
gray = cr_gt(data1,sigma,res,max_dat)
# [150:230,150:230]
fig = px.imshow(gray,color_continuous_scale='spectral',width=400, height=400)
fig.update_layout(coloraxis_showscale=False)
# fig.update_layout( title=Patinet11)
cols_p2[0].plotly_chart(fig)

idx=np.where(sample_id==Patinet12+'_'+Array1+'_'+Cycle1+'_'+Exposure_time1)[0][0]
data1=kinome_data[sample_id1[idx]].astype('float')
gray = cr_gt(data1,sigma,res,max_dat)
# [150:230,150:230]
fig = px.imshow(gray,color_continuous_scale='spectral',width=400, height=400)
fig.update_layout(coloraxis_showscale=False)
# fig.update_layout( title=Patinet12)
cols_p2[1].plotly_chart(fig)

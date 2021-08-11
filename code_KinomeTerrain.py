# -*- coding: utf-8 -*-
"""
@author: Ehsan Saghapour
"""
import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
import plotly.express as px
import networkx as nx
import plotly.figure_factory as ff
from scipy.spatial.distance import pdist, squareform
import seaborn as sns
from matplotlib import pyplot as plt

data= pd.read_csv('Kinome_SIGNAL_(MinBknd).csv',sep=',',header=None)
data1= pd.read_csv('log2_kinome.csv',sep=',',header=None)

ID_peptide=data1.iloc[8:,0]

peptic_name=data.iloc[11:,0]
patient_id=np.unique(data.iloc[5,7:])
Array=np.unique(data.iloc[2,7:])
Cycle=np.unique(data.iloc[3,7:])
Exposure_time=np.unique(data.iloc[4,7:])
Drug=[ 'CT', 'Brigatinib', 'Sitravatinib','Neratinib']
Cycle_plt=np.array(data.iloc[3,7:75])
Exposure_time_plt=np.array(data.iloc[4,7:75])
Cycle_signal=pd.DataFrame(np.zeros((68,3)))
Cycle_signal.columns=['Cycle','Signal','Exposure_Time']
Cycle_signal['Cycle']=Cycle_plt
Cycle_signal['Exposure_Time']=Exposure_time_plt


sample_id=data.iloc[5,7:]+'_'+data.iloc[2,7:] +'_'+data.iloc[3,7:].astype('str')+'_'+data.iloc[4,7:].astype('str')

kinome_data=pd.DataFrame(data.iloc[11:,7:])
kinome_data.columns=sample_id
kinome_data.index=peptic_name



A1= np.zeros((21, 194 ,68))

itr=0
for i in range(0+7,(68*4*21)+7,68*4):
    A1[itr,:,:]= data.iloc[11:,i:i+68]
    itr=itr+1

A11=np.median(A1,axis=0)


A2= np.zeros((21, 194 ,68))

itr=0
for i in range(68+7,(68*4*21)+7,68*4):
    A2[itr,:,:]= data.iloc[11:,i:i+68]
    itr=itr+1

A3= np.zeros((21, 194 ,68))

itr=0
for i in range(2*68+7,(68*4*21)+7,68*4):
    A3[itr,:,:]= data.iloc[11:,i:i+68]
    itr=itr+1
    
A4= np.zeros((21, 194 ,68))

itr=0
for i in range(3*68+7,(68*4*21)+7,68*4):
    A4[itr,:,:]= data.iloc[11:,i:i+68]
    itr=itr+1


st.title('KinomeTerrain (Version 1)')
st.text('By E.Saghapour, J.Anderson, K.Lee, C.Willey(Leader).')

# st.header('Plot Cycle/Exposure_Time to Signal for each Peptide in different arrays and Tumors')
st.header('Kinetic Phosphorylation Over Time')


col1,col2,col3= st.beta_columns(3)
with col1:
    peptide1 = st.selectbox(
                "ID_peptide", ID_peptide.to_list(),key=11)
               
with col2:
    Array1 = st.selectbox(
                "Array", Array,key=11 )
with col3:
    Tumor1= st.selectbox(
                "Tumor", patient_id,key=11 )
    
# with col4:
#     plt_cycle_exposure= st.selectbox(
#                 "Cycle/Exposure", ['Cycle','Exposure_Time'])

plt_cycle_exposure ="Cycle"

agree = st.checkbox('Exposure_Time (Default is Cycle)')
if agree: 
    plt_cycle_exposure='Exposure_Time'
    

idx=np.where(data.iloc[5,7:]==Tumor1)[0][0]//(68*4)
idx_peptide=np.where(peptic_name==peptide1)[0][0]
print(idx_peptide)
print(idx)
btn1= st.checkbox('Run')

# col1,col2= st.beta_columns(2)

def cr_two_fig(Cycle_signal):
    Cycle_signal=Cycle_signal.astype('float')
    Cycle_signal['Cycle']=Cycle_signal['Cycle'].astype('int')
    Cycle_signal1=Cycle_signal[Cycle_signal['Cycle']!=154]
    
    print(Cycle_signal1)
    fig = px.line(Cycle_signal1, x="Exposure_Time", y="Signal", color="Cycle")
    fig.update_traces(mode='markers+lines')  

    Cycle_signal1=Cycle_signal[Cycle_signal['Cycle']==154]
    
    print(Cycle_signal1)
    fig1 = px.line(Cycle_signal1, x="Exposure_Time", y="Signal", color="Cycle")
    fig1.update_traces(mode='markers+lines') 
    # st.plotly_chart(fig)
    # st.plotly_chart(fig1)
    fig.update_layout(
        autosize=False,
        width=350,
        height=400,

        )
    fig1.update_layout(
        autosize=False,
        width=350,
        height=400,

        )
    col1,col2= st.beta_columns(2)
    with col1:
        st.plotly_chart(fig)
    with col2:
        st.plotly_chart(fig1)
    
if btn1:
    if Array1=='A1':
        Cycle_signal['Signal']=A1[idx,idx_peptide,:]
        mm=np.mean(A1[idx,idx_peptide,-8:-1])
        mm1=np.mean(A1[idx,idx_peptide,0:5])

        if plt_cycle_exposure=='Cycle':
            fig = px.scatter(Cycle_signal, x="Cycle", y="Signal", color='Exposure_Time')
            fig.add_trace(go.Scatter(x=[154, 32],y=[mm,mm1], mode="lines", name='Slope'))
            st.plotly_chart(fig)
            # fig1 = px.scatter(Cycle_signal, x="Cycle", y="Signal", color='Exposure_Time')
            # fig1.add_trace(go.Scatter(x=[154, 32],y=[mm,mm1], mode="lines", name='Slope'))
            # fig.update_traces(mode='markers+lines')
        else:
            cr_two_fig(Cycle_signal)
    
    elif Array1=='A2':
        Cycle_signal['Signal']=A2[idx,idx_peptide,:]
        mm=np.mean(A2[idx,idx_peptide,-8:-1])
        mm1=np.mean(A2[idx,idx_peptide,0:5])
        if plt_cycle_exposure=='Cycle':
            fig = px.scatter(Cycle_signal, x="Cycle", y="Signal", color='Exposure_Time')
            fig.add_trace(go.Scatter(x=[154, 32],y=[mm,mm1], mode="lines", name='Slope'))
        else:
            cr_two_fig(Cycle_signal)
            
        st.plotly_chart(fig)
    elif Array1=='A2':
        Cycle_signal['Signal']=A3[idx,idx_peptide,:]
        mm=np.mean(A3[idx,idx_peptide,-8:-1])
        mm1=np.mean(A3[idx,idx_peptide,0:5])
        if plt_cycle_exposure=='Cycle':
            fig = px.scatter(Cycle_signal, x="Cycle", y="Signal", color='Exposure_Time')
            fig.add_trace(go.Scatter(x=[154, 32],y=[mm,mm1], mode="lines", name='Slope'))
            st.plotly_chart(fig)
        else:
            cr_two_fig(Cycle_signal)  
            
        
        
    
    else:
        Cycle_signal['Signal']=A4[idx,idx_peptide,:]
        # Cycle_signal=Cycle_signal.sort_values(by=['Cycle'],ascending=False)
        mm=np.mean(A4[idx,idx_peptide,-8:-1])
        mm1=np.mean(A4[idx,idx_peptide,0:5])
        if plt_cycle_exposure=='Cycle':
            fig = px.scatter(Cycle_signal, x="Cycle", y="Signal", color='Exposure_Time')
            fig.add_trace(go.Scatter(x=[154, 32],y=[mm,mm1], mode="lines", name='Slope'))
            st.plotly_chart(fig)
        else:
            cr_two_fig(Cycle_signal)  


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

# st.title('KinomeTerrain')
st.text('******************************************************************************************')


st.header('Drug Response Terrain Patterning')

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

btn3= st.checkbox('Run',key=3)
if btn3:
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
    
    
                # fig.update_layout(coloraxis_showscale=False)
                fig.update_layout( title=Array[itr]+ ' ( ' +Drug[itr]+' )')
                cols[i].plotly_chart(fig)
                itr=itr+1
              

# st.plotly_chart(node_trace)
st.text('******************************************************************************************')

st.header('Cross-tumor Comparative Terrain Modeling')

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
    
# cols1,cols2= st.beta_columns(2)

# Patinet11 = st.multiselect(patient_id)

Patinet11 = st.multiselect("Tumors",patient_id)
print(len(Patinet11))
btn4= st.checkbox('Run',key= 4)
if btn4:
    cols = st.beta_columns(len(Patinet11)) 
    
    sample_id1=sample_id.to_list()
    itr=0
    max_dat=[]
    for j in range(len(Patinet11)):
                print(Patinet11[j]+'_'+Array1+'_'+Cycle1+'_'+Exposure_time1)
                print(sample_id)
                idx=np.where(sample_id==Patinet11[j]+'_'+Array1+'_'+Cycle1+'_'+Exposure_time1)
                print(idx[0][0])
                data1=kinome_data[Patinet11[j]+'_'+Array1+'_'+Cycle1+'_'+Exposure_time1].astype('float')
                max_dat.append(np.max(data1))
                
                
    print(max(max_dat))
    max_dat=max(max_dat)
    
    
    itr=0
    if len(Patinet11)%2 ==1: k=np.round(len(Patinet11)/2)+1 
    else:  k=np.round(len(Patinet11)/2)
    for j in range(int(k)):
            cols = st.beta_columns(2) 
            # re.columns=['MSE','SSIM']
            for i in range(2):
                if len(Patinet11) >= itr:
                    print(Patinet11[itr]+'_'+Array1+'_'+Cycle1+'_'+Exposure_time1)
                    # idx=np.where(sample_id=Patinet11[itr]+'_'+Array1+'_'+Cycle1+'_'+Exposure_time1)[0][0]
                    data1=kinome_data[Patinet11[itr]+'_'+Array1+'_'+Cycle1+'_'+Exposure_time1].astype('float')
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
        
        
                    # fig.update_layout(coloraxis_showscale=False)
                    fig.update_layout( title=Patinet11[itr])
                    cols[i].plotly_chart(fig)
                    itr=itr+1



def cr_heat_denogram(data_array,labels,labels_col,title):      #This fucntion had gotten from https://plotly.com/python/dendrogram/
    # Initialize figure by creating upper dendrogram
    fig = ff.create_dendrogram(data_array, orientation='bottom', labels=labels)
    for i in range(len(fig['data'])):
        fig['data'][i]['yaxis'] = 'y2'
    
    # for data in fig['data']:
    #     fig.add_trace(data)
    
    # Create Side Dendrogram
    dendro_side = ff.create_dendrogram(data_array.T, orientation='right', labels=labels_col)
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'
    
    # Add Side Dendrogram Data to Figure
    for data in dendro_side['data']:
        fig.add_trace(data)
    
    # Create Heatmap
    # dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
    # dendro_leaves = list(map(int, dendro_leaves))
    # data_dist = pdist(data_array.T)
    # heat_data = squareform(data_dist)
    # heat_data = heat_data[dendro_leaves,:]
    # heat_data = heat_data[:,dendro_leaves]
    
    heatmap = [
        go.Heatmap(
            x = labels,
            y = labels_col,
            z = data_array.T,
            colorscale = 'Blues'
        )
    ]
    
    heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
    heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']
    
    # Add Heatmap Data to Figure
    for data in heatmap:
        fig.add_trace(data)
    
    # Edit Layout
    fig.update_layout({'width':800, 'height':800,
                              'showlegend':False, 'hovermode': 'closest',
                              })
    # Edit xaxis
    fig.update_layout(xaxis={'domain': [.15, 1],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'ticks':""})
    # Edit xaxis2
    fig.update_layout(xaxis2={'domain': [0, .15],
                                        'mirror': False,
                                        'showgrid': False,
                                        'showline': False,
                                        'zeroline': False,
                                        'showticklabels': False,
                                        'ticks':""})
    
    # Edit yaxis
    fig.update_layout(yaxis={'domain': [0, .85],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'showticklabels': False,
                                      'ticks': ""
                            })
    # Edit yaxis2
    fig.update_layout(yaxis2={'domain':[.825, .975],
                                        'mirror': False,
                                        'showgrid': False,
                                        'showline': False,
                                        'zeroline': False,
                                        'showticklabels': False,
                                        'ticks':""})
    fig.update_layout(title=title)
    #               yaxis={'labels_col.to_list'})
    # Plot!
    return fig

def cr_heat_denogram1(data_array,labels,title):      #This fucntion had gotten from https://plotly.com/python/dendrogram/
    # Initialize figure by creating upper dendrogram
    fig = ff.create_dendrogram(data_array, orientation='bottom', labels=labels)
    for i in range(len(fig['data'])):
        fig['data'][i]['yaxis'] = 'y2'
    
    # Create Side Dendrogram
    dendro_side = ff.create_dendrogram(data_array, orientation='right')
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'
    
    # Add Side Dendrogram Data to Figure
    for data in dendro_side['data']:
        fig.add_trace(data)
    
    # Create Heatmap
    dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
    dendro_leaves = list(map(int, dendro_leaves))
    data_dist = pdist(data_array)
    heat_data = squareform(data_dist)
    heat_data = heat_data[dendro_leaves,:]
    heat_data = heat_data[:,dendro_leaves]
    
    heatmap = [
        go.Heatmap(
            x = dendro_leaves,
            y = dendro_leaves,
            z = heat_data,
            colorscale = 'Blues'
        )
    ]
    
    heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
    heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

    # Add Heatmap Data to Figure
    for data in heatmap:
        fig.add_trace(data)
    
    # Edit Layout
    fig.update_layout({'width':800, 'height':800,
                              'showlegend':False, 'hovermode': 'closest',
                              })
    # Edit xaxis
    fig.update_layout(xaxis={'domain': [.15, 1],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'ticks':""})
    # Edit xaxis2
    fig.update_layout(xaxis2={'domain': [0, .15],
                                        'mirror': False,
                                        'showgrid': False,
                                        'showline': False,
                                        'zeroline': False,
                                        'showticklabels': False,
                                        'ticks':""})
    
    # Edit yaxis
    fig.update_layout(yaxis={'domain': [0, .85],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'showticklabels': False,
                                      'ticks': ""
                            })
    # Edit yaxis2
    fig.update_layout(yaxis2={'domain':[.825, .975],
                                        'mirror': False,
                                        'showgrid': False,
                                        'showline': False,
                                        'zeroline': False,
                                        'showticklabels': False,
                                        'ticks':""})
    fig.update_layout(title=title)
    #               yaxis={'labels_col.to_list'})
    # Plot!
    return fig
# st.plotly_chart(node_trace)
st.text('******************************************************************************************')

st.header('Plot a Dendrogram with a Heatmap for Arrays')

col1,col2= st.beta_columns(2)
with col1:
    Array1 = st.selectbox(
                "Array", ['A1','A2','A3','A4','All_Array'],key=2
            )
# get data
data1= pd.read_csv('log2_kinome.csv',sep=',',header=None)
labels_col=data1.iloc[1,7:91]+'_'+data1.iloc[2,7:91]
labels=data1.iloc[8:,0]
labels=labels.to_list()
st.set_option('deprecation.showPyplotGlobalUse', False)
btn5= st.checkbox('Run',key=5)
if btn5:
    if Array1 != 'All_Array':
        idx=np.where(data1.iloc[2,0:91]==Array1)[0]
        data_array=data1.iloc[8:,idx]
        labels_col=data1.iloc[1,idx]+'_'+data1.iloc[2,idx]
        data_array=data_array.values
    
        fig1=cr_heat_denogram1(data_array,labels,'Peptide vs Peptide')
        st.plotly_chart(fig1)
        fig2=cr_heat_denogram1(data_array.T,labels_col.to_list(),'Tumor vs Tumor')
        st.plotly_chart(fig2)
        fig3=cr_heat_denogram(data_array,labels,labels_col.to_list(),'Peptide vs Tumor')
        st.plotly_chart(fig3)
        data_array=pd.DataFrame(data_array).astype('float')
        data_array.index=labels
        data_array.columns=labels_col.to_list()
        
        # sns.clustermap(data_array,standard_scale=1, metric="correlation",figsize=(10, 30))
        # sns.clustermap(data_array, standard_scale=1)
        sns.clustermap(data_array, standard_scale=1,figsize=(20, 30))
        st.pyplot()
        sns.clustermap(data_array.T, standard_scale=1,figsize=(30, 20))
        st.pyplot()
    else:
        
        data_array=data1.iloc[8:,7:91]
        data_array=data_array.values
        labels_col=data1.iloc[1,7:91]+'_'+data1.iloc[2,7:91]
        fig1=cr_heat_denogram1(data_array,labels,'Peptide vs Peptide')
        st.plotly_chart(fig1)
        fig2=cr_heat_denogram1(data_array.T,labels_col.to_list(),'Tumor vs Tumor')
        st.plotly_chart(fig2)
        fig3=cr_heat_denogram(data_array,labels,labels_col.to_list(),'Peptide vs Tumor')
        st.plotly_chart(fig3)
        data_array=pd.DataFrame(data_array).astype('float')
        data_array.index=labels
        data_array.columns=labels_col.to_list()
        
        # sns.clustermap(data_array,standard_scale=1, metric="correlation",figsize=(10, 30))
        # sns.clustermap(data_array, standard_scale=1)
        sns.clustermap(data_array, standard_scale=1,figsize=(20, 30))
        st.pyplot()
        sns.clustermap(data_array.T, standard_scale=1,figsize=(30, 20))
        st.pyplot()

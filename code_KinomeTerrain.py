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
from sklearn.linear_model import LinearRegression

def app():
    
    data= pd.read_csv('Kinome_SIGNAL_(MinBknd).csv',sep=',',header=None,encoding= 'unicode_escape')
    data1= pd.read_csv('log2_kinome.csv',sep=',',header=None,encoding= 'unicode_escape')
    st.sidebar.text("-------------------------------------")

    # ID_peptide=data1.iloc[8:,0]
    
    peptic_name=data.iloc[11:,0]
    patient_id=np.unique(data.iloc[5,7:])
    Array=np.unique(data.iloc[2,7:])
    Cycle=np.unique(data.iloc[3,7:])
    # Exposure_time=np.unique(data.iloc[4,7:])
    Drug=[ 'CT', 'Brigatinib', 'Sitravatinib','Neratinib']
    Cycle_plt=np.array(data.iloc[3,7:75])
    Exposure_time_plt=np.array(data.iloc[4,7:75])
    Cycle_signal=pd.DataFrame(np.zeros((68,3)))
    Cycle_signal.columns=['Cycle','Signal','Exposure_Time']
    Cycle_signal['Cycle']=Cycle_plt
    Cycle_signal['Exposure_Time']=Exposure_time_plt
    
    
    sample_id=data.iloc[5,7:]+'_'+data.iloc[2,7:] +'_'+data.iloc[3,7:].astype('str')+'_'+data.iloc[4,7:].astype('str')
    
    kinome_data=pd.DataFrame(data.iloc[11:,7:]).copy()    
    
    from sklearn.preprocessing import StandardScaler
    
    scaler = StandardScaler()
    


    for i in range(0,(68*4*21),68*4):
        data_t=kinome_data.iloc[:,i:i+272].values.T
        scaler.fit(data_t)
        kinome_data.iloc[:,i:i+272]=scaler.transform(data_t).T
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
    
    
    st.title('KinomeTerrain (Version 3)')
    st.text('By E.Saghapour, J.Anderson, J.Chen(Leader), C.Willey(Leader).')
    
    
    
    
    
    # st.header('Plot Cycle/Exposure_Time to Signal for each Peptide in different arrays and Tumors')
    st.header('Kinetic Phosphorylation Over Time')
    with st.expander('Importance peptide for each Tumor'):
        data2=data.copy()
        col1,col2,col3,col4= st.columns(4)
        with col2:
            Tmr = st.selectbox(
                        "Tumor", patient_id,key=8)
            idx=np.where(data2.iloc[5,7:]==Tmr)[0]
        with col1:
            num_input= st.number_input('The number of importance peptide',value=5,min_value =2)
            
            mat_tumor=data2.iloc[11:,idx[:-1]]
            print('*********')
            print(mat_tumor)
            max_mat_tmr= np.max(mat_tumor,axis=1)
            max_mat_tmr=pd.DataFrame(max_mat_tmr,columns=['max_all_col'])
            max_mat_tmr.index=peptic_name.to_list()
            max_mat_tmr.index.name = 'Importance_peptide'
            max_mat_tmr['max_all_col']=max_mat_tmr['max_all_col'].astype('int')
            max_mat_tmr=max_mat_tmr.sort_values(by=['max_all_col'], ascending=False)
            print(num_input)
            st.write(max_mat_tmr.index[0:int(num_input)])
    col1,col2,col3,col4= st.columns(4)
    with col1:
        peptide1 = st.selectbox(
                    "ID_peptide", peptic_name.to_list(),key=11)
                   
    with col2:
        Array1 = st.selectbox(
                    "Array", Array,key=11 )
    with col3:
        Tumor1= st.selectbox(
                    "Tumor", patient_id,key=11 )
        
    with col4:
        plt_cycle_exposure= st.selectbox(
                    "Cycle/Exposure", ['Cycle','Exposure_Time'])
    
    # plt_cycle_exposure ="Cycle"
    
    # agree = st.checkbox('Exposure_Time (Default is Cycle)')
    # if agree: 
    #     plt_cycle_exposure='Exposure_Time'
        
    
    idx=np.where(data.iloc[5,7:]==Tumor1)[0][0]//(68*4)
    idx_peptide=np.where(peptic_name==peptide1)[0][0]
    # print(idx_peptide)
    # print(idx)
    btn1= st.checkbox('Run')
    
    # col1,col2= st.columns(2)
    
    def cr_two_fig(Cycle_signal):
        Cycle_signal=Cycle_signal.astype('float')
        Cycle_signal['Cycle']=Cycle_signal['Cycle'].astype('int')
        Cycle_signal1=Cycle_signal[Cycle_signal['Cycle']!=154]
        
        #print(Cycle_signal1)
        fig = px.line(Cycle_signal1, x="Exposure_Time", y="Signal", color="Cycle")
        fig.update_traces(mode='markers+lines')  
    
        Cycle_signal1=Cycle_signal[Cycle_signal['Cycle']==154]
        
        X = Cycle_signal1.Exposure_Time.values[:-1].reshape(-1, 1)
        model = LinearRegression()
        model.fit(X, Cycle_signal1.Signal[:-1])
        
        x_range = np.linspace(5, 205, 50)
        y_range = model.predict(x_range.reshape(-1, 1))
        fig1 = px.line(Cycle_signal1, x="Exposure_Time", y="Signal", color="Cycle")
        fig1.update_traces(mode='markers+lines') 
 
        fig1.add_traces(go.Scatter(x=x_range, y=y_range, name='Expected line',
                                   line = dict(shape = 'linear', width= 2, dash = 'dash')))
                        
        fig.update_layout(
            autosize=False,
            width=600,
            height=400,
    
            )
        fig1.update_layout(
            autosize=False,
            width=600,
            height=400,
    
            )
        col1,col2= st.columns(2)
        with col1:
            st.plotly_chart(fig)
        with col2:
            st.plotly_chart(fig1)
            
    def plt_overtime(A):
        Cycle_signal['Signal']=A[idx,idx_peptide,:]
        mm=np.mean(A[idx,idx_peptide,-8:-1])
        mm1=np.mean(A[idx,idx_peptide,0:5])
        if plt_cycle_exposure=='Cycle':
            fig = px.scatter(Cycle_signal, x="Cycle", y="Signal", color='Exposure_Time')
            fig.add_trace(go.Scatter(x=[154, 32],y=[mm,mm1], mode="lines", name='Slope'))
            st.plotly_chart(fig)
        else:
            cr_two_fig(Cycle_signal)
    if btn1:
        if Array1=='A1':
           plt_overtime(A1)
        
        elif Array1=='A2':
            plt_overtime(A2)
                
        elif Array1=='A3':
            plt_overtime(A3) 
    
        else:
            plt_overtime(A4) 
    
    
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
    
    
    
    
    
    
    # sigma=.04
    
    # sample_id=sample_id.to_list()
    def graph(name):
            name_graph=['Kamada-Kawai','circular',
                         'random','fruchterman_reingold','bipartite','spring','shell'] 
            idx=np.where(np.array(name_graph)==name)[0]
            #print(name)
            #print(idx)
            if idx==0:
                pos= nx.kamada_kawai_layout(G)
            elif idx==1:
                pos= nx.circular_layout(G)
    
            elif idx==2:
                pos= nx.random_layout(G)
            elif idx==3:
                pos= nx.fruchterman_reingold_layout(G)
            elif idx==4:
                pos= nx.bipartite_layout(G,G)
            elif idx==5:
                pos= nx.spring_layout(G)
            elif idx==6:
                pos= nx.shell_layout(G)
            pos=nx.rescale_layout_dict(pos)
            return (pos)
                
    
    
    
    def cr_gt(data1,sigma,res,max_dat,pos):  # Create for Gene terrain 
        
        x_ = np.linspace(-1, 1, res)
        y_ = np.linspace(-1, 1, res)
        X, Y = np.meshgrid(x_, y_)
        
        gaussian=np.zeros((res,res))
    
        pi=3.14
        itr=0
        for node in G.nodes():
            x, y = pos[node]
            amp1=data1[itr]
            # amp1=np.exp(data1[itr]/max_dat)

            itr=itr+1
            gaussian1 =(1/np.sqrt(2*pi*sigma))*np.exp(-(((X-x)/sigma)**2+((Y-y)/sigma)**2)) 
            gaussian =amp1* gaussian1 + gaussian
            # gaussian = np.flipud(gaussian)
            # gray=(gaussian-np.min(gaussian))/(np.max(gaussian)-np.min(gaussian))
        
        return gaussian
    
    # st.title('KinomeTerrain')
    st.text('***********************************************************************************************************************************')
    
    
    st.header('Drug Response Terrain Patterning')
    
    col1,col2,col3,col4= st.columns(4)
    with col1:
        Patinet1 = st.selectbox(
                    "Tumor", patient_id
                )
    
    with col2:
    
        Cycle1 = st.selectbox(
                    "Cycle", Cycle
                )
    if Cycle1 != str(154): 
        with col3:
        
            Exposure_time1 = st.selectbox(
                        "Exposure_time", ['20','50','100']
                    )
    else:
        with col3:
            Exposure_time1 = st.selectbox(
                        "Exposure_time", ['10','20','50','100','200']
                    )
    
    with col4:
        M_layout = st.selectbox(
                    "Position nodes(layout)", ['Kamada-Kawai','circular',
                         'random','bipartite','spring']
                )
    pos=graph(M_layout)
    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = pos[node]
        node_x.append(200*x+200)
        node_y.append(200*y+200)
    nme=[]
    for node in G.nodes():
        nme.append(node)
    nme=[]
    for node in G.nodes():
        nme.append(node)
    
    
    st.sidebar.header('Properties of Terrain')
    sigma=st.sidebar.slider('Sigma',0.0,0.1,.04,.01)
    btn= st.sidebar.checkbox('W/Peptide names')
    # res=st.sidebar.slider('Resulotion of image',100, 1000, 400,100)
    res=400
    btn3= st.checkbox('Run',key=3)
    if btn3:
        cols = st.columns(4) 
        
        sample_id1=sample_id.to_list()

        itr=0
        for j in range(2):
                cols = st.columns(2) 
                # re.columns=['MSE','SSIM']
                for i in range(2):
                    print(Patinet1+'_'+Array[j]+'_'+Cycle1+'_'+Exposure_time1)
                    idx=np.where(sample_id==Patinet1+'_'+Array[itr]+'_'+Cycle1+'_'+Exposure_time1)[0][0]
                    data1=kinome_data[sample_id1[idx]].astype('float')
                    gray = cr_gt(data1,sigma,res,1,pos)
                    # [150:230,150:230]
                    # fig = px.imshow(gray,color_continuous_scale='spectral',width=500, height=500)
                    fig = px.imshow(gray,color_continuous_scale='jet',range_color=[-3,3],width=500, height=500)
                    if btn:
                        
                        fig.add_trace(go.Scatter(x=node_x, y=node_y,
                                  text=nme,
                                  textposition='top center',
                        mode='markers',
                        name='Point'
                                  ))
        
        
                    # fig.update_layout(coloraxis_showscale=False)
                    fig.update_layout( title=Array[itr]+ ' ( ' +Drug[itr]+' )')
                    cols[i].plotly_chart(fig)
                    itr=itr+1
                  
    
    # st.plotly_chart(node_trace)
    st.text('***********************************************************************************************************************************')
    
    st.header('Cross-tumor Comparative Terrain Modeling')
    
    print(Cycle)
    col1,col2,col3,col4= st.columns(4)
    with col1:
        Array1 = st.selectbox(
                    "Array", Array,key=1
                )
    
    with col2:
        
        Cycle1 = st.selectbox(
                    "Cycle", Cycle ,key=1
                )
        
    if Cycle1 != str(154): 
        with col3:
        
            Exposure_time1 = st.selectbox(
                        "Exposure_time", ['20','50','100'], key=1
                    )
    else:
        with col3:
            Exposure_time1 = st.selectbox(
                        "Exposure_time", ['10','20','50','100','200'],key=1
                    )
    with col4:
        M_layout = st.selectbox(
                    "Position nodes(layout)", ['Kamada-Kawai','circular',
                         'random','bipartite','spring'],key=1
                )    
    pos=graph(M_layout)
    
    
    # cols1,cols2= st.columns(2)
    
    # Patinet11 = st.multiselect(patient_id)
    
    Patinet11 = st.multiselect("Tumors",patient_id)
    
    btn4= st.checkbox('Run',key= 4)
    if btn4:
    
        if len(Patinet11)>1:
            cols = st.columns(len(Patinet11)) 
        
            sample_id1=sample_id.to_list()
            
            itr=0
            if len(Patinet11)%2 ==1: k=np.round(len(Patinet11)/2)+1 
            else:  k=np.round(len(Patinet11)/2)
            for j in range(int(k)):
                    cols = st.columns(2) 
                    # re.columns=['MSE','SSIM']
                    for i in range(2):
                        if len(Patinet11) >= itr:
                            print(Patinet11[itr]+'_'+Array1+'_'+Cycle1+'_'+Exposure_time1)
                            # idx=np.where(sample_id=Patinet11[itr]+'_'+Array1+'_'+Cycle1+'_'+Exposure_time1)[0][0]
                            data1=kinome_data[Patinet11[itr]+'_'+Array1+'_'+Cycle1+'_'+Exposure_time1].astype('float')
                            gray = cr_gt(data1,sigma,res,1,pos)
                            # [150:230,150:230]
                            # fig = px.imshow(gray,color_continuous_scale='spectral',width=500, height=500)
                            fig = px.imshow(gray,color_continuous_scale='jet',range_color=[-3,3],width=500, height=500)
                            if btn:
                                
                                fig.add_trace(go.Scatter(x=node_x, y=node_y,
                                          text=nme,
                                          textposition='top center',
                                mode='markers',
                                name='Point'
                                          ))
                
                
                            # fig.update_layout(coloraxis_showscale=False)
                            fig.update_layout( title=Patinet11[itr])
                            cols[i].plotly_chart(fig)
                            itr=itr+1
        else:
            st.warning('Please select at least two tumors')
    
    
    def cr_heat_denogram(df,title):      #This fucntion had gotten from https://plotly.com/python/dendrogram/
            # Initialize figure by creating upper dendrogram
              dataHeat_arr= df     
              dataHeat_arr_t= np.transpose(dataHeat_arr)
              fig = ff.create_dendrogram(dataHeat_arr_t, orientation='bottom', labels=df.columns  )


              for i in range(len(fig['data'])):
                  fig['data'][i]['yaxis'] = 'y2'

              # Create Side Dendrogram

              # dendro_side = ff.create_dendrogram(dataHeat_arr_t, orientation='right' ,labels=name_molec[:100])
              dendro_side = ff.create_dendrogram(dataHeat_arr, orientation='right', labels=df.index)
              for i in range(len(dendro_side['data'])):
                  dendro_side['data'][i]['xaxis'] = 'x2'


              # Add Side Dendrogram Data to Figure
              for data in dendro_side['data']:
                  fig.add_trace(data)

              heatmap = [
                  go.Heatmap(
                      x = df.index ,
                      y =df.columns ,
                      z = dataHeat_arr,
                      colorscale = 'Blues',

                      colorbar = dict(
                              title="Scale",
                              thicknessmode="pixels",
                              thickness=20,
                              yanchor="top",
                              tickfont=dict(size=10),
                              x=0,
                              y=1,
                              len=.1,

                              )
                  )
              ]

              heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
              heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

              # Add Heatmap Data to Figure
              for data in heatmap:
                  fig.add_trace(data)

              fig['layout']['yaxis']['ticktext'] = np.asarray(df.index)
              fig['layout']['yaxis']['tickvals'] = np.asarray(dendro_side['layout']['yaxis']['tickvals'])

              # Edit Layout
              fig.update_layout({'width':800, 'height':1000,
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
                                                'showticklabels': True,
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
              fig.update_layout(
                  yaxis={'side': 'right'} ,
              )
              fig.update_layout(plot_bgcolor='white')
              return fig
    
    def cr_heat_denogram1(data_array,labels,title):      #This fucntion is avaibale through https://plotly.com/python/dendrogram/
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
        fig.update_layout(
                  yaxis={'side': 'right'} ,
              )
        return fig
    # st.plotly_chart(node_trace)
    st.text('***********************************************************************************************************************************')
    
    st.header('Plot a Dendrogram with a Heatmap for Arrays')
    
    col1,col2= st.columns(2)
    with col1:
        Array1 = st.selectbox(
                    "Array", ['A1','A2','A3','A4','All_Array'],key=2
                )
    # get data
    data1= pd.read_csv('log2_kinome.csv',sep=',',header=None,encoding= 'unicode_escape')
    labels_col=data1.iloc[1,7:91]+'_'+data1.iloc[2,7:91]
    labels_tumor=[]
    itr=0
    for name in data1.iloc[1,7:91]:
         idx_tmr=np.where(data.iloc[1,:]==name)[0][0]
         labels_tumor.append(data.iloc[5,idx_tmr]+'_'+data1.iloc[2,7+itr])
         itr=itr+1
     
    labels=data1.iloc[8:,0]
    labels=labels.to_list()
    st.set_option('deprecation.showPyplotGlobalUse', False)
    btn5= st.checkbox('Run',key=5)
    if btn5:
        if Array1 != 'All_Array':
            idx=np.where(data1.iloc[2,0:91]==Array1)[0]
            data_array=data1.iloc[8:,idx]
            labels_col=data1.iloc[1,idx]+'_'+data1.iloc[2,idx]
            
            labels_tumor=[]
            itr=0
            for name in labels_col:
                 idx_tmr=np.where(data.iloc[1,:]==name[0:-3])[0][0]
                 labels_tumor.append(data.iloc[5,idx_tmr]+'_'+Array1)
                 itr=itr+1
            
            
            data_array=data_array.values.astype('float')
            
            fig1=cr_heat_denogram1(data_array,labels,'Peptide vs Peptide')
            st.plotly_chart(fig1)
            fig2=cr_heat_denogram1(data_array.T,labels_tumor,'Tumor vs Tumor')
            st.plotly_chart(fig2)
            
            data_array=pd.DataFrame(data_array).astype('float')
            data_array.index=labels
            data_array.columns=labels_tumor
            
            fig3=cr_heat_denogram(data_array,'Peptide vs Tumor')
            st.plotly_chart(fig3)
            
            # sns.clustermap(data_array,standard_scale=1, metric="correlation",figsize=(10, 30))
            # sns.clustermap(data_array, standard_scale=1)
            sns.clustermap(data_array, standard_scale=1,figsize=(20, 30))
            st.pyplot()
            sns.clustermap(data_array.T, standard_scale=1,figsize=(30, 20))
            st.pyplot()
        else:
            
            data_array=data1.iloc[8:,7:91].astype('float')
            # data_array=data_array.values
            labels_col=data1.iloc[1,7:91]+'_'+data1.iloc[2,7:91]
            fig1=cr_heat_denogram1(data_array,labels,'Peptide vs Peptide')
            st.plotly_chart(fig1)
            fig2=cr_heat_denogram1(data_array.T,labels_tumor,'Tumor vs Tumor')
            st.plotly_chart(fig2)
            
            data_array=pd.DataFrame(data_array).astype('float')
            data_array.index=labels
            data_array.columns=labels_tumor
            fig3=cr_heat_denogram(data_array,'Peptide vs Tumor')
            st.plotly_chart(fig3)
            # sns.clustermap(data_array,standard_scale=1, metric="correlation",figsize=(10, 30))
            # sns.clustermap(data_array, standard_scale=1)
            sns.clustermap(data_array, standard_scale=1,figsize=(20, 30))
            st.pyplot()
            sns.clustermap(data_array.T, standard_scale=1,figsize=(30, 20))
            st.pyplot()

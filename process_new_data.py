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
import base64
from sklearn.linear_model import LinearRegression

def app():

    st.title('KinomeTerrain (Version 2.1)')
    st.text('By E.Saghapour, J.Anderson, Jake Chen(Leader), C.Willey(Leader).')
    
    st.header('Analysis your data !')
    
    
    st.sidebar.text("-------------------------------------")
    st.sidebar.write("If you do not have any data, downlaod a sample. Then, upload it into the app")
    btn5=st.sidebar.button('Click Here ! ')
    
    def get_table_download_link_csv(df):

            # df = df)
            # df=df.drop(['index'],axis=1)
            csv = df.to_csv(header=False, index=False, index_label=None)
            b64 = base64.b64encode(csv.encode()).decode()  # some strings <-> bytes conversions necessary here
            href = f'<a href="data:file/csv;base64,{b64}" download="Toy_example.csv" >Download ! </a>'
            st.sidebar.markdown(href, unsafe_allow_html=True)
            
    if btn5 :
        Gene_expr=pd.read_csv('Toy_example.csv',sep=',',header=None )
        get_table_download_link_csv(Gene_expr)
    st.sidebar.text("-------------------------------------")

    def file_selector():
        file = st.file_uploader("Choose a file", type="csv")
        if file is not None:
          # Gene_expr=pd.read_csv(file,sep='\t')
          Gene_expr=pd.read_csv(file,sep=',',header=None )
          return Gene_expr
    
    data=file_selector()
    def Info_data(Gene_expr):
        st.header( "Information about File:")
        st.write( "The Number of Rows and Columns are:" , np.shape(Gene_expr)  , ", respectively. " )
        # st.write( "The number of Genes is:" , np.shape(Gene_expr)  )
        st.write('The 20 First Rows of your data:', Gene_expr.head(20))
    
    if data is not None:
        Info_data(data)
        
        # data= pd.read_csv('Toy_example.csv',sep=',',header=None)
        
        peptic_name=data.iloc[4:,0]
        # patient_id=np.unique(data.iloc[5,7:])
        Array=np.unique(data.iloc[0,1:])
        Cycle=np.unique(data.iloc[1,1:])
        # Exposure_time=np.unique(data.iloc[2,1:])
        Drug=[ 'CT', 'Brigatinib', 'Sitravatinib','Neratinib']
        Cycle_plt=np.array(data.iloc[1,1:69])
        Exposure_time_plt=np.array(data.iloc[2,1:69])
        Cycle_signal=pd.DataFrame(np.zeros((68,3)))
        Cycle_signal.columns=['Cycle','Signal','Exposure_Time']
        Cycle_signal['Cycle']=Cycle_plt
        Cycle_signal['Exposure_Time']=Exposure_time_plt
        
        
        sample_id=data.iloc[0,1:]+'_'+data.iloc[1,1:].astype('str')+'_'+data.iloc[2,1:].astype('str')
        
        kinome_data=pd.DataFrame(data.iloc[4:,1:])
        kinome_data.columns=sample_id
        kinome_data.index=peptic_name
        
        
        
        
        # st.header('Plot Cycle/Exposure_Time to Signal for each Peptide in different arrays and Tumors')
        st.header('Kinetic Phosphorylation Over Time')
        
        
        col1,col2,col3,= st.beta_columns(3)
        with col1:
            peptide1 = st.selectbox(
                        "ID_peptide", peptic_name.to_list(),key=11)
                       
        with col2:
            Array1 = st.selectbox(
                        "Array", Array,key=11 )
            
        with col3:
            plt_cycle_exposure= st.selectbox(
                        "Cycle/Exposure", ['Cycle','Exposure_Time'])
        
        
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
            
            X = Cycle_signal1.Exposure_Time.values[:-1].reshape(-1, 1)
            model = LinearRegression()
            model.fit(X, Cycle_signal1.Signal[:-1])
            
            x_range = np.linspace(5, 205, 50)
            y_range = model.predict(x_range.reshape(-1, 1))
            fig1 = px.line(Cycle_signal1, x="Exposure_Time", y="Signal", color="Cycle")
            fig1.update_traces(mode='markers+lines') 
     
            fig1.add_traces(go.Scatter(x=x_range, y=y_range, name='Expected linear',
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
            col1,col2= st.beta_columns(2)
            with col1:
                st.plotly_chart(fig)
            with col2:
                st.plotly_chart(fig1)
        idx_peptide=np.where(peptic_name==peptide1)[0][0]
        print(idx_peptide)
        A1= np.array(data.iloc[4:,1:69])
        A2= np.array(data.iloc[4:,69:2*68+1])
        A3= np.array(data.iloc[4:,2*68+1:3*68+1])
        A4= np.array(data.iloc[4:,3*68+1:4*68+1]) 
        if btn1:
            if Array1=='A1':
                Cycle_signal['Signal']=A1[idx_peptide,:]
                mm=np.mean(A1[idx_peptide,-8:-1].astype(int))
                mm1=np.mean(A1[idx_peptide,0:5].astype(int))
        
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
                Cycle_signal['Signal']=A2[idx_peptide,:].astype(int)
                mm=np.mean(A2[idx_peptide,-8:-1].astype(int))
                mm1=np.mean(A2[idx_peptide,0:5].astype(int))
                if plt_cycle_exposure=='Cycle':
                    fig = px.scatter(Cycle_signal, x="Cycle", y="Signal", color='Exposure_Time')
                    fig.add_trace(go.Scatter(x=[154, 32],y=[mm,mm1], mode="lines", name='Slope'))
                else:
                    cr_two_fig(Cycle_signal)
                    
            elif Array1=='A3':
                Cycle_signal['Signal']=A3[idx_peptide,:].astype(int)
                mm=np.mean(A3[idx_peptide,-8:-1].astype(int))
                mm1=np.mean(A3[idx_peptide,0:5].astype(int))
                if plt_cycle_exposure=='Cycle':
                    fig = px.scatter(Cycle_signal, x="Cycle", y="Signal", color='Exposure_Time')
                    fig.add_trace(go.Scatter(x=[154, 32],y=[mm,mm1], mode="lines", name='Slope'))
                    st.plotly_chart(fig)
                else:
                    cr_two_fig(Cycle_signal)  
        
            else:
                Cycle_signal['Signal']=A4[idx_peptide,:].astype(int)
                mm=np.mean(A4[idx_peptide,-8:-1].astype(int))
                mm1=np.mean(A4[idx_peptide,0:5].astype(int))
                if plt_cycle_exposure=='Cycle':
                    fig = px.scatter(Cycle_signal, x="Cycle", y="Signal", color='Exposure_Time')
                    fig.add_trace(go.Scatter(x=[154, 32],y=[mm,mm1], mode="lines", name='Slope'))
                    st.plotly_chart(fig)
                else:
                    cr_two_fig(Cycle_signal)  
        
        data1= pd.read_csv('Kinome_SIGNAL_(MinBknd).csv',sep=',',header=None)
        A1= np.zeros((21, 194 ,68))
        
        itr=0
        for i in range(0+7,(68*4*21)+7,68*4):
            A1[itr,:,:]= data1.iloc[11:,i:i+68]
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
        # # pos= nx.kamada_kawai_layout(G)
        # # pos=nx.random_layout(G)
        # pos=nx.spring_layout(G)
        # pos=nx.rescale_layout_dict(pos)
        # G.nodes()
        
        
        
        
        # sigma=.04
        
        # sample_id=sample_id.to_list()
        def graph(name):
                name_graph=['Kamada-Kawai','circular',
                             'random','fruchterman_reingold','bipartite','spring','shell'] 
                idx=np.where(np.array(name_graph)==name)[0]
                print(name)
                print(idx)
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
                amp1=np.exp(data1[itr]/max_dat)
                itr=itr+1
                gaussian1 =(1/np.sqrt(2*pi*sigma))*np.exp(-(((X-x)/sigma)**2+((Y-y)/sigma)**2)) 
                gaussian =amp1* gaussian1 + gaussian
                # gaussian = np.flipud(gaussian)
                gray=(gaussian-np.min(gaussian))/(np.max(gaussian)-np.min(gaussian))
            
            return 1-gray
        
        # st.title('KinomeTerrain')
        st.text('***********************************************************************************************************************************')
        
        
        st.header('Drug Response Terrain Patterning')
        
        col1,col2,col3,col4= st.beta_columns(4)
        
        with col1:
        
            Cycle1 = st.selectbox(
                        "Cycle", Cycle
                    )
        if Cycle1 != str(154): 
            with col2:
            
                Exposure_time1 = st.selectbox(
                            "Exposure_time", ['20','50','100']
                        )
        else:
            with col2:
                Exposure_time1 = st.selectbox(
                            "Exposure_time", ['10','20','50','100','200']
                        )
        
        with col3:
            M_layout = st.selectbox(
                        "Position nodes(layout)", ['Kamada-Kawai','circular',
                             'random','bipartite','spring']
                    )
        pos=graph(M_layout)
        
        
        
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
                        print(Array[j]+'_'+Cycle1+'_'+Exposure_time1)
                        idx=np.where(sample_id==Array[j]+'_'+Cycle1+'_'+Exposure_time1)[0][0]
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
                        print(Array[j]+'_'+Cycle1+'_'+Exposure_time1)
                        idx=np.where(sample_id==Array[itr]+'_'+Cycle1+'_'+Exposure_time1)[0][0]
                        data1=kinome_data[sample_id1[idx]].astype('float')
                        gray = cr_gt(data1,sigma,res,max_dat,pos)
                        # [150:230,150:230]
                        fig = px.imshow(gray,color_continuous_scale='spectral',width=500, height=500)
                        if btn:
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
                            nme=[]
                            for node in G.nodes():
                                nme.append(node)
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

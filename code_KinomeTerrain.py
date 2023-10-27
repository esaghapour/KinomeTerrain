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
            
    def plt_overtime(A,idx):
        Cycle_signal['Signal']=A[idx,idx_peptide,:]
        mm=np.mean(A[idx,idx_peptide,-8:-1])
        mm1=np.mean(A[idx,idx_peptide,0:5])
        if plt_cycle_exposure=='Cycle':
            fig = px.scatter(Cycle_signal, x="Cycle", y="Signal", color='Exposure_Time')
            fig.add_trace(go.Scatter(x=[154, 32],y=[mm,mm1], mode="lines", name='Slope'))
            st.plotly_chart(fig)
        else:
            cr_two_fig(Cycle_signal)
    
    
    
    data= pd.read_csv('Kinome_SIGNAL_(MinBknd).csv',sep=',',header=None,encoding= 'unicode_escape')
    data1= pd.read_csv('log2_kinome.csv',sep=',',header=None,encoding= 'unicode_escape')
    Layout= pd.read_csv('Layout_Kinome.txt',sep=',')
    Layout=Layout.set_index('name')
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
    
    # from sklearn.preprocessing import StandardScaler
    
    # scaler = StandardScaler()
    
    kinome_data=kinome_data.astype('float')

    for i in range(0,(68*4*21),68*4):
        # data_t=kinome_data.iloc[:,i:i+272].values.T
        # scaler.fit(data_t)
        kinome_data.iloc[:,i:i+272]=(kinome_data.iloc[:,i:i+272].values)-np.mean(kinome_data.iloc[:,i:i+272].values)/np.std(kinome_data.iloc[:,i:i+272].values)
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
    
    
    st.title('KinomeTerrain (Version 4)')
    st.text('By E.Saghapour, J.Anderson, J.Chen(Leader), C.Willey(Leader).')
    
    st.text('')
    st.text('')
    tab1, tab2, tab3 = st.tabs(["Kinetic Phosphorylation Over Time",
                                      "Kinease Response in Terrain", "HeatMap"])

    
    
    

    
    

    with tab1:
        st.header('Kinetic Phosphorylation Over Time')
        with st.expander('Importance peptide for each Tumor'):
            
            data2=data.copy()
            data2=pd.DataFrame(data.iloc[1:,6:].values)
            data2.iloc[9:,0]=data.iloc[10:,0]
            data2=data2.set_index(0)
        
            data2 = data2.rename(index={0: 'index'})
            
            col1,col2,col3,col4= st.columns(4)
            with col2:
                Tmr = st.selectbox(
                            "Tumor", patient_id,key=8)
                idx=np.where(data2.loc['Sample name',:]==Tmr)[0]
            with col1:
                num_input= st.number_input('The number of importance peptide',value=5,min_value =2)
                
                
                mat_tumor=data2.iloc[10:,idx].astype('float')
                
              
                print('*********')
                print(mat_tumor)
                max_mat_tmr= np.max(mat_tumor,axis=1)
                max_mat_tmr=pd.DataFrame(max_mat_tmr,columns=['Max_value'])
                max_mat_tmr=max_mat_tmr.sort_values(by=['Max_value'], ascending=False)
                max_mat_tmr.index.name = 'Importance peptide'
                max_mat_tmr=max_mat_tmr.astype(int)
                print(num_input)
                st.write(max_mat_tmr.iloc[0:int(num_input),:])
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
        
        
        if btn1:
            if Array1=='A1':
               plt_overtime(A1,idx)
            
            elif Array1=='A2':
                plt_overtime(A2,idx)
                    
            elif Array1=='A3':
                plt_overtime(A3,idx)
        
            else:
                plt_overtime(A4,idx)
        
        
        
        ################## Show individual peptides based on exposure time or cycle  ###########
        st.subheader('Show Individual Peptides Based On Exposure Time or Cycle')
        data2=data.copy()
        data2=pd.DataFrame(data.iloc[1:,6:].values)
        data2.iloc[9:,0]=data.iloc[10:,0]
        data2=data2.set_index(0)
    
        data2 = data2.rename(index={0: 'index'})
        col1,col2,col3,col4= st.columns(4)
        with col1:
            Tumor1= st.selectbox(
                        "Tumor", patient_id,key=21 )
                       
        with col2:
            Array1 = st.selectbox(
                        "Array", Array,key=21 )
        with col3:
            cycle_exposure= st.selectbox(
                        "Cycle/Exposure", ['Cycle','Exposure_Time'],key=21)
            
        with col4:
            if cycle_exposure=='Cycle':
                Cycle1 = st.selectbox(
                            "Cycle", np.unique(Cycle_plt),key=21
                        )
            else: 
                Exposure_time1 = st.selectbox(
                            "Exposure_time", ['10','20','50','100','200'],key=21
                        )
            
        peptide1 = st.multiselect(
                    "ID_peptide", peptic_name.to_list(),key=21)    
            
            
        
        btn21= st.checkbox('Run',key=21)
    
        
        
        if btn21:
            
            idx=np.where(data2.loc['Sample name',:]==Tumor1)[0]
            data3=data2.iloc[:,idx]
            idx=np.where(data3.loc['Array',:]==Array1)[0]
            data3=data3.iloc[:,idx]
            
            if cycle_exposure=='Cycle':
                idx=np.where(data3.loc['Cycle',:]==Cycle1)[0]
                
                data3=data3.iloc[:,idx]
                expr_time=data3.loc['Exposure time',:]
                data3=data3.loc[peptide1,:].astype('float')
                
                data3=pd.DataFrame(data3.values,columns=expr_time.values)
                
                data3.index=peptide1
                print(data3.T)
                fig = px.line(data3.T, markers=True)
                fig.update_layout(xaxis_title="Exposure time",yaxis_title="Signal value",legend_title="Peptide name")
    
                st.plotly_chart(fig)
                
            else:
                idx=np.where(data3.loc['Exposure time',:]==Exposure_time1)[0]
                data3=data3.iloc[:,idx]
                val_cycle=data3.loc['Cycle',:]
                data3=data3.loc[peptide1,:].astype('float')
                
                data3=pd.DataFrame(data3.values,columns=val_cycle.values)
                
                data3.index=peptide1
                print(data3.T)
                fig = px.line(data3.T, markers=True)
                fig.update_layout(xaxis_title="Cycle",yaxis_title="Signal value",legend_title="Peptide name")
                st.plotly_chart(fig)
        
    
                
    
    
    
    def cr_gt(data1,Layout,sigma,res,name):  # Create for Gene terrain 
        
        x_ = np.linspace(-1, 1, res)
        y_ = np.linspace(-1, 1, res)
        X, Y = np.meshgrid(x_, y_)
        
        gaussian=np.zeros((res,res))
    
        pi=3.14

        print(data1)
        print(Layout)
        hlp_data=Layout.merge(data1,how='inner',left_index=True,right_index=True)
        print(hlp_data)
        # res=200


              # print(name)
        gaussian=np.zeros((res,res))
  
        for i in range(np.shape(hlp_data)[0]):
            # print(i)
            x=hlp_data['X'][i]
            y=hlp_data['Y'][i]
            amp1=hlp_data[name][i]
            
            gaussian1 =(1/np.sqrt(2*pi*sigma))*np.exp(-(((X-x)/sigma)**2+((Y-y)/sigma)**2)) 
            gaussian =amp1* gaussian1 + gaussian
        
        
        return gaussian
    

    with tab2:
        st.header('Kinease Response')
        
        col1,col2,col3= st.columns(3)
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
        
        
        
        st.sidebar.header('Properties of Terrain')
        sigma=st.sidebar.slider('Sigma',0.0,0.1,.04,.01)
        btn= st.sidebar.checkbox('W/Peptide names')
        raw_mark = dict(
                color='black',    # Changed marker color to look cleaner with a boundary
                size=5,
                symbol='circle',
                opacity=0.4,      # Brought down the opacity and boundaries often allow for this  # We give the pre-defined dictionary to this attribute.
        )
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
                        print('***********')
                        print(Patinet1+'_'+Array[itr]+'_'+Cycle1+'_'+Exposure_time1)
                        print('***********')
                        idx=np.where(sample_id==Patinet1+'_'+Array[itr]+'_'+Cycle1+'_'+Exposure_time1)[0][0]
                        data1=kinome_data[sample_id1[idx]].astype('float')
                        print(data1)
                        gray = cr_gt(data1,Layout,sigma,res,Patinet1+'_'+Array[itr]+'_'+Cycle1+'_'+Exposure_time1)
                        # [150:230,150:230]
                        # fig = px.imshow(gray,color_continuous_scale='spectral',width=500, height=500)
                        import plotly.graph_objects as go

                        fig = go.Figure(data=go.Heatmap(z=gray, colorscale='Jet', zmin=-100, zmax=100))
                        fig.update_layout(width=500, height=500)
                        if btn:
                            
                            fig.add_trace(go.Scatter(x=200*Layout['X']+200, y=200*Layout['Y']+200,
                                      text=Layout.index,
                                      textposition='top center',
                            mode='markers',
                            name='Point',
                            textfont=dict(
                           
                            size=10,
                            color="blue"
                        ),
                            marker=raw_mark,
    
    
                            
                                      ))
            
            
                        # fig.update_layout(coloraxis_showscale=False)
                        fig.update_layout( title=Array[itr]+ ' ( ' +Drug[itr]+' )')
                        cols[i].plotly_chart(fig)
                        itr=itr+1
                      
        
        # st.plotly_chart(node_trace)
        
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

        # pos=graph(M_layout)
        
        
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
                                gray = cr_gt(data1,Layout,sigma,res,Patinet11[itr]+'_'+Array1+'_'+Cycle1+'_'+Exposure_time1)
                                # [150:230,150:230]
                                # fig = px.imshow(gray,color_continuous_scale='spectral',width=500, height=500)
                                fig = px.imshow(gray,color_continuous_scale='jet',range_color=[-100,100],width=500, height=500)
                                if btn:
                                    
                                    fig.add_trace(go.Scatter(x=200*Layout['X']+200, y=200*Layout['Y']+200,
                                              text=Layout.index,
                                              textposition='top center',
                                    mode='markers',
                                    name='Point',
                                    textfont=dict(
                                   
                                    size=10,
                                    color="blue"
                                ),
                                    marker=raw_mark,
            
            
                                    
                                              ))
                    
                    
                                # fig.update_layout(coloraxis_showscale=False)
                                fig.update_layout( title=Patinet11[itr])
                                cols[i].plotly_chart(fig)
                                itr=itr+1
            else:
                st.warning('Please select at least two tumors')
    
    
    def cr_heat_denogram(df,title, width=800, height= 1000,Color='jet'):      #This fucntion had gotten from https://plotly.com/python/dendrogram/
      # Initialize figure by creating upper dendrogram 10*20
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

      fig_leaves = fig['layout']['xaxis']['ticktext']
      #fig_leaves = list(map(int, fig_leaves))

      dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
      #dendro_leaves = list(map(int, dendro_leaves))


      dataHeat_arr = dataHeat_arr.loc[dendro_leaves,:]
      dataHeat_arr = dataHeat_arr.loc[:,fig_leaves]


      heatmap = [
          go.Heatmap(
              x = fig_leaves ,
              y =dendro_leaves ,
              z = dataHeat_arr,
              colorscale = Color,
              colorbar = dict(
                      title="Scale",
                      thicknessmode="pixels",
                      thickness=20,
                      yanchor="top",
                      tickfont=dict(size=10),
                      x=1.2,
                      y=1,
                      len=.2,

                      )
          )
      ]

      heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
      heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

      # Add Heatmap Data to Figure
      for data in heatmap:
          fig.add_trace(data)

      fig['layout']['yaxis']['ticktext'] = np.asarray(dendro_side['layout']['yaxis']['ticktext'])
      fig['layout']['yaxis']['tickvals'] = np.asarray(dendro_side['layout']['yaxis']['tickvals'])

      # Edit Layout
      fig.update_layout({'width':width, 'height':height,
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

      fig.update_layout(
          yaxis={'side': 'right'} ,
      )
      fig.update_layout(plot_bgcolor='white')
      return fig
    

    # st.plotly_chart(node_trace)
    with tab3:
        st.header('Plot a Dendrogram with a Heatmap for Arrays')
        
        col1,col2,col3= st.columns(3)
        with col1:
            Array1 = st.selectbox(
                        "Array", ['A1','A2','A3','A4','All_Array'],key=5
                    )
            
        with col2:
            Color = st.selectbox(
                        "Color", ['blues','reds','jet','Picnic'],key=5
                    )
        with col3:
            metric = st.selectbox(
                        "Metric (works for peptides or tummors)", ['Euclidean','Correlation'],key=5
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
                
                if metric=='Euclidean':
                    
                    data_dist = pdist(data_array)
                    heat_data = squareform(data_dist)
                else:
                    heat_data = np.corrcoef(data_array)
               
                heat_data=pd.DataFrame(heat_data)
                heat_data.index=labels
                heat_data.columns=labels
    
                fig1=cr_heat_denogram(heat_data,'Peptide vs Peptide',800,800,Color)
                st.plotly_chart(fig1)
    
    
                if metric=='Distance':
                    
                    data_dist = pdist(data_array.T)
                    heat_data = squareform(data_dist)
                else:
                    heat_data = np.corrcoef(data_array.T)
    
                heat_data=pd.DataFrame(heat_data)
                heat_data.index=labels_tumor
                heat_data.columns=labels_tumor
    
                fig2=cr_heat_denogram(heat_data,'Tumor vs Tumor',800, 800,Color)
                st.plotly_chart(fig2)
                
                data_array=pd.DataFrame(data_array).astype('float')
                data_array.index=labels
                data_array.columns=labels_tumor
                
    
                fig3=cr_heat_denogram(data_array,'Peptide vs Tumor',800, 1000,Color)
                st.plotly_chart(fig3)

                sns.clustermap(data_array, standard_scale=1,figsize=(20, 30))

            else:
                
                data_array=data1.iloc[8:,7:91].astype('float')
                # data_array=data_array.values
                labels_col=data1.iloc[1,7:91]+'_'+data1.iloc[2,7:91]
                
                data_dist = pdist(data_array)
                heat_data = squareform(data_dist)
                heat_data=pd.DataFrame(heat_data)
                heat_data.index=labels
                heat_data.columns=labels
    
                fig1=cr_heat_denogram(heat_data,'Peptide vs Peptide',800,800,Color)
                st.plotly_chart(fig1)
    
                data_dist = pdist(data_array.T)
                heat_data = squareform(data_dist.T)
                heat_data=pd.DataFrame(heat_data)
                heat_data.index=labels_tumor
                heat_data.columns=labels_tumor
    
                fig2=cr_heat_denogram(heat_data,'Tumor vs Tumor',800,800,Color)
                st.plotly_chart(fig2)
                
                data_array=pd.DataFrame(data_array).astype('float')
                data_array.index=labels
                data_array.columns=labels_tumor
                
                
                fig3=cr_heat_denogram(data_array,'Peptide vs Tumor',800,1000,Color)
                st.plotly_chart(fig3)

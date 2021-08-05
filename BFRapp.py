#!/usr/bin/env python
# coding: utf-8

# In[2]:


#packages
import streamlit as st
import base64
from Bio import SeqIO
from Bio.SeqUtils import GC, molecular_weight as mw,MeltingTemp as mt
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import io

# In[ ]:
st.set_option('deprecation.showPyplotGlobalUse', False)

# app logo
_,col2,_ =st.beta_columns([0.6,8,0.6])
with col2:
    logo=open("logo/logo.png","rb")
    l=logo.read()
    dat=base64.b64encode(l).decode("UTF-8")
    logo.close()
    st.markdown(f'<img src="data:image/png;base64,{dat}" width="100%" alt="logo png">',unsafe_allow_html=True)



# Menu        
Menu=["About","DNA Lab","Proteins Lab"]    
menu=st.sidebar.selectbox("Menu",Menu)
    
# DNA Lab
if menu=="About":
    # About this Application
    st.markdown("***")
    st.subheader("About This App")
    st.markdown("")
    st.write("**Bioinformatics Files Reader** or **BFR** is a streamlit web application that can read and process common bioinformatics files. BFR supports **Fasta** files, **GenBank** files and **GenPept** files.")
    st.write("This application is made with:  **Streamlit**, **Biopython**, **Numpy**, **Matplotlib** and **Pandas**.")
    # About me
    st.markdown("***")
    st.subheader("About The App Developer")
    st.markdown("")
    #my name
    st.write("**Full Name: **")
    st.write("*Yassine OUCHEN*")
    #contact
    st.write("**Contact: **")
    
    
    # my gmail
    
    col1,col2=st.beta_columns([1.2,20])
    with col1:
        logo=open("logo/gmail.png","rb") 
        l=logo.read()
        dat=base64.b64encode(l).decode("UTF-8")
        logo.close()
        st.markdown(f'<img src="data:image/png;base64,{dat}" width="100%" alt="gmail png">',unsafe_allow_html=True)
    with col2:
        st.markdown('<a href="yassine.ouchen2021@gmail.com">yassine.ouchen2021@gmail.com</a>',unsafe_allow_html=True)
        
    # my linkedin
    
    col1,col2=st.beta_columns([1.2,20])
    with col1:
        logo=open("logo/link.png","rb") 
        l=logo.read()
        dat=base64.b64encode(l).decode("UTF-8")
        logo.close()
        st.markdown(f'<img src="data:image/png;base64,{dat}" width="100%" alt="link png">',unsafe_allow_html=True)
    with col2:
        st.write("https://www.linkedin.com/in/yassine-ouchen-9772a81aa")
    
    # my github
    
    col1,col2=st.beta_columns([1.2,20])
    with col1:
        logo=open("logo/github.png","rb") 
        l=logo.read()
        dat=base64.b64encode(l).decode("UTF-8")
        logo.close()
        st.markdown(f'<img src="data:image/png;base64,{dat}" width="100%" alt="github png">',unsafe_allow_html=True)
    with col2:
        st.write("https://github.com/ouchen-bio")
        
    # references
    st.sidebar.markdown("***")
    col1,col2=st.sidebar.beta_columns([3,20])
    with col1:
        logo=open("logo/ref.png","rb") 
        l=logo.read()
        dat=base64.b64encode(l).decode("UTF-8")
        logo.close()
        st.markdown(f'<img src="data:image/png;base64,{dat}" width="100%" alt="ref png">',unsafe_allow_html=True)
    with col2:
        st.write("**References**")
        
    st.markdown("")
    st.sidebar.write("**C.N. Pace, F. Vajdos, L. Fee, G. Grimsley, T. Grey**. *How to measure and predict the molar absorption coefficient of a protein? Protein Sci, 4,(1995), pp. 2411-2423.*")
    st.markdown("")
    st.sidebar.write("**Guruprasad K, Reddy BV, Pandit MW (1990)**. *Correlation between stability of a protein and its dipeptide composition: a novel approach for predicting in vivo stability of a protein from its primary sequence. Protein Eng. 4 (2): 155–61.*")
    st.markdown("")
    st.sidebar.write("**Kyte, J. and Doolittle, R. 1982.** *A simple method for displaying the hydropathic character of a protein. J. Mol. Biol. 157: 105-132.*")
    st.markdown("")
    st.sidebar.write("**Lobry, J. R., and C. Gautier, (1994)**.*Hydrophobicity, expressivity and aromaticity are the major trends of amino acid usage in 999 Escherichia coli chromosome encoded genes. Nucleic Acids Research 22: 3174-3180.*")
    
    
                    
# DNA Lab        
elif menu=="DNA Lab":
    st.markdown("***")
    st.sidebar.markdown("***")
    col1,col2,col3=st.sidebar.beta_columns([2,14,2])
    with col2:
        st.write("**Welcome To DNA Labortory !**")
    #DNA Gif Photo
    gifbar1=open("logo/dna3.gif","rb")
    cont=gifbar1.read()
    data1=base64.b64encode(cont).decode("UTF-8")
    gifbar1.close()
    st.sidebar.markdown(f'<img src="data:image/gif;base64,{data1}" width="100%" alt="dna gif">',unsafe_allow_html=True)
    st.sidebar.markdown("***")
    st.sidebar.write("**Note:** Please make sure to upload a **Fasta** file or a **GenBank** file that containts **DNA** data. ")
    # DNA Lab Content
    st.header("File Analysis")
    file_uploaded=st.file_uploader("Upload Your File",type=["fasta","gb"])
    # CONVERTING FASTA FORMAT INTO A READABLE STREAMLIT FILE FARMAT
    if file_uploaded is not None:
        st.markdown("***")
        byte_str=file_uploaded.read()
        text_obj=byte_str.decode("UTF-8")
        
        if file_uploaded.name[-5:]=="fasta":
            seq_object=SeqIO.read(io.StringIO(text_obj),"fasta")
        else:
            seq_object=SeqIO.read(io.StringIO(text_obj),"gb")
        
        # file_description
        file_description=seq_object.description
        st.subheader("***1- File Description***")
        st.write(file_description)
        st.markdown("")
        
        # DNA Sequence
        st.subheader("***2- DNA Sequence***")
        DNA_seq=seq_object.seq
        generate_seq=st.button("Generate DNA Sequence")
        if generate_seq is True:
            st.write(DNA_seq)
        st.markdown("")
        
        # statistics
        st.subheader("***3- Statistics***")
        st.markdown("")
        
        # DNA Length
        DNA_length=len(DNA_seq)
        st.markdown("")
        
        # total number of each nucleotide
        Adenine=DNA_seq.count("A")
        Thymine=DNA_seq.count("T")
        Cytosine=DNA_seq.count("C")
        Guanine=DNA_seq.count("G")
        # table_count
        def DNA_Numbers(sequence):
            d=dict([("Adenine",Adenine),
             ("Thymine",Thymine),
             ("Cytosine",Cytosine),
             ("Guanine",Guanine)
            ])
            return d
        X1=DNA_Numbers(DNA_seq)
        count_table=pd.DataFrame.from_dict(X1,orient="index")
        count_table=count_table.rename({0:"Count"},axis="columns")
        count_table.reset_index(inplace=True)
        count_table=count_table.rename(columns={"index":"Nucleotide"})
        
        
        # Nucleotide Frequency
        Adeninep=(DNA_seq.count("A")/DNA_length)*100
        Thyminep=(DNA_seq.count("T")/DNA_length)*100
        Cytosinep=(DNA_seq.count("C")/DNA_length)*100
        Guaninep=(DNA_seq.count("G")/DNA_length)*100
        
        # table_percentage
        def DNA_Numbers(sequence):
            d=dict([("Adenine",Adeninep),
             ("Thymine",Thyminep),
             ("Cytosine",Cytosinep),
             ("Guanine",Guaninep)
            ])
            return d
        X2=DNA_Numbers(DNA_seq)
        percent_table=pd.DataFrame.from_dict(X2,orient="index")
        percent_table=percent_table.rename({0:"Percentage %"},axis="columns")
        percent_table.reset_index(inplace=True)
        percent_table=percent_table.rename(columns={"index":"Nucleotide"})
        
        
        # gc content
        GC=GC(DNA_seq)
        
        st.write("**3.1- DNA Length**")
        
        st.write(DNA_length,"Nucleotides")
        st.markdown("")
        
        st.write("**3.2- Nucleotide Frequency**")
        st.markdown("")
        
        RD=st.radio("Select",["Numbers","Percentages"])
        if RD=="Numbers":
            st.write(count_table)
            st.markdown("")
            Dotplot1=st.radio("Dotplot",["No","Yes"])
            if Dotplot1=="Yes":
              # dotplot number
                dotn=plt.bar(X1.keys(),X1.values())
                st.pyplot()
            
        else:
            st.write(percent_table)
            st.markdown("")
            Dotplot2=st.radio("Dotplot",["No","Yes"])
            if Dotplot2=="Yes":   
                # dotplot percent
                dotP=plt.bar(X2.keys(),X2.values())
                st.pyplot()
         
        st.markdown("")    
        # gc content
        st.write("**3.3- GC Content**")
        st.write(GC,"%")
        st.markdown("")
        
        # Physico-Chemical 
        st.write("***4- Physico-Chemical Properties***")
        st.markdown("")
        
        # molecular Weight
        st.write("**4.1- Molecular Weight**")
        molec_weight=mw(DNA_seq,"DNA")
        st.write("MW≈",molec_weight,"Da")
        st.markdown("")
        
        # melting temperature
        st.write("**4.2- Melting Temperature**")
        MT_radio=st.radio("Select",["Standard","Other Factors"])
        
        # standard melting temperature
        def m_t_s(seq):
            d=dict([("Tm_GC","%0.2f" % mt.Tm_GC(seq)),
             ("Tm_NN","%0.2f" % mt.Tm_NN(seq))
            ])
            return d
        T1=m_t_s(DNA_seq)
        MTS_table=pd.DataFrame.from_dict(T1,orient="index")
        MTS_table=MTS_table.rename({0:"Value (°C)"},axis="columns")
        MTS_table.reset_index(inplace=True)
        MTS_table=MTS_table.rename(columns={"index":"Method"})
       
  
        if MT_radio=="Standard":
            st.write(MTS_table)
        else:
            st.markdown("")
            # other Factors
            na=st.slider("Na (mMol)",0,100,50)
            mg=st.slider("Mg (mMol)",0.0,3.0,1.5)
            tris=st.slider("Tris (mMol)",0,20,10)
            dntps=st.slider("dNTPs (mMol)",0.0,1.0,0.5)
            salt=st.selectbox("salt correlation",[0,1,2,3,4,6,7])
            Tm_GCp="%0.2f" % mt.Tm_GC(DNA_seq,Na=na,Mg=mg,Tris=tris,dNTPs=dntps,saltcorr=salt)
            Tm_NNp="%0.2f" % mt.Tm_NN(DNA_seq,Na=na,Mg=mg,Tris=tris,dNTPs=dntps,saltcorr=salt)
        
            def m_t_p(seq):
                d=dict([("Tm_GC",Tm_GCp),
                 ("Tm_NN",Tm_NNp)
                ])
                return d
            T2=m_t_p(DNA_seq)
            MTP_table=pd.DataFrame.from_dict(T2,orient="index")
            MTP_table=MTP_table.rename({0:"Value (°C)"},axis="columns")
            MTP_table.reset_index(inplace=True)
            MTP_table=MTP_table.rename(columns={"index":"Method"})
            st.markdown("")
            st.write(MTP_table)

# protein Lab        
elif menu=="Proteins Lab":
    st.markdown("***")
    st.sidebar.markdown("***")
    col1,col2,col3=st.sidebar.beta_columns([2,20,2])
    with col2:
        st.write("**Welcome To Proteins Labortory !**")
    #Protein Gif Photo
    gifbar1=open("logo/protein.gif","rb")
    cont=gifbar1.read()
    data1=base64.b64encode(cont).decode("UTF-8")
    gifbar1.close()
    st.sidebar.markdown(f'<img src="data:image/gif;base64,{data1}" width="100%" alt="protein gif">',unsafe_allow_html=True)
    st.sidebar.markdown("***")
    st.sidebar.write("**Note:** Please make sure to upload a **Fasta** file, **GenBank** file or **GenPept** file that containts **protein** data. ")
    # Protein Lab Content
    st.header("File Analysis")
    file_name=st.file_uploader("Upload Your File",type=["fasta","gp","gb"])
    
    # CONVERTING FASTA and GenPept FORMAT INTO A READABLE STREAMLIT FILE FARMAT
    if file_name is not None:
        st.markdown("***")
        byte_str=file_name.read()
        text_obj=byte_str.decode("UTF-8")
        
        if file_name.name[-5:]=="fasta":
            seq_object_protein=SeqIO.read(io.StringIO(text_obj),"fasta")
        else:
            # gb file=gp file
            seq_object_protein=SeqIO.read(io.StringIO(text_obj),"gb")
        
        # file_description
        file_description=seq_object_protein.description
        st.subheader("1- File Description")
        st.write(file_description)
        st.markdown("")
        
        # generate Amino Acids sequences
        st.subheader("2- Protein Sequence")
        AA_seq=seq_object_protein.seq
        generate_seq=st.button("Generate Protein Sequence")
        if generate_seq is True:
            st.write(AA_seq)
        st.markdown("")
        
        # apply ProtParam module(see pckgs)
        # turn AA_seq into string
        protein_seq=str(AA_seq)
        protein_analysed=ProteinAnalysis(protein_seq)
        
        # statistics
        st.subheader("3- Statistics")
        st.markdown("")
        
        st.write("**3.1- Protein Lenght**")
        P_lenght=len(AA_seq)
        st.write(P_lenght,"aa")
        st.markdown("")
        
        st.write("**3.2- Amino Acids Count**")
        AA_count=protein_analysed.count_amino_acids() # return a dict of number of each AA
        
        AA_Dotplot=st.radio("Dotplot",["No","Yes"])
        
        # amino acids count Table
        AA_table=pd.DataFrame.from_dict(AA_count,orient="index")
        AA_table=AA_table.rename({0:"Count"},axis="columns")
        AA_table.reset_index(inplace=True)
        AA_table=AA_table.rename(columns={"index":"Amino Acid"})
        st.write(AA_table)
        # dotplot
        if AA_Dotplot=="Yes":
            st.markdown("")
            dotP=plt.bar(AA_count.keys(),AA_count.values())
            st.pyplot()
        st.markdown("")
        
        st.subheader("4- Physico-Chemical Properties")
        st.markdown("")
        
        # molecular weight
        st.write("**4.1- Molecular Weight**")
        p_MW=("%0.2f" % protein_analysed.molecular_weight())
        st.write(float(p_MW),"Da")
        st.markdown("")
        
        # Aromaticity according to Lobry,1994
        st.write("**4.2- Aromaticity** (Lobry,1994)")
        p_Arm=("%0.2f" % protein_analysed.aromaticity())
        st.write(float(p_Arm))
        st.markdown("")
        
        # instability index
        st.write("**4.3- Instability Index** (Guruprasad et al.1990)")
        p_ii=("%0.2f" % protein_analysed.instability_index())
        st.write(float(p_ii))
        st.markdown("")
        
        # Isoelectric point 
        st.write("**4.3- Isoelectric Point**")
        p_pi=("%0.2f" % protein_analysed.isoelectric_point())
        st.write(float(p_pi))
        st.markdown("")
        
        # Protein charge at given pH
        st.write("**4.3- Charge**")
        pH=st.slider("Select pH value",0.0,14.0,float(p_pi))
        p_cpH=("%0.2f" % protein_analysed.charge_at_pH(pH))
        st.write(float(p_cpH))
        st.markdown("")
        
        # Gravy according to Kyte and Doolittle
        st.write("**4.5- GRAVY **(Kyte and Doolittle)")
        p_gv=("%0.2f" % protein_analysed.gravy())
        st.write(float(p_gv))
        st.markdown("")
        
        # estimate the molar extinction coefficient
        st.write("**4.6- Molar Extinction Coefficient**")
        p_mec=protein_analysed.molar_extinction_coefficient()
        
        def MEC(cof):
            d=dict([("with reduced cysteines",cof[0]),
             ("with disulfid bridges",cof[1])
            ])
            return d
        C=MEC(p_mec)
        MEC_table=pd.DataFrame.from_dict(C,orient="index")
        MEC_table=MEC_table.rename({0:"Molar_Extinction_Coefficient (/M.cm)"},axis="columns")
        MEC_table.reset_index(inplace=True)
        MEC_table=MEC_table.rename(columns={"index":"Variable"})
        st.write(MEC_table)
        st.markdown("")
        
        # Secondary Stracture Fractions
        st.write("**4.7- Secondary Stracture Fraction**")
        sec_struc=protein_analysed.secondary_structure_fraction()
    
        def SSF(fact):
            d=dict([("Helix","%0.2f" % fact[0]),
             ("Turn","%0.2f" % fact[1]),
             ("Sheet","%0.2f" % fact[2])
            ])
            return d
        S=SSF(sec_struc)
        SSF_table=pd.DataFrame.from_dict(S,orient="index")
        SSF_table=SSF_table.rename({0:"Fraction"},axis="columns")
        SSF_table.reset_index(inplace=True)
        SSF_table=SSF_table.rename(columns={"index":"Part "})
        st.write(SSF_table)
       
        
        
               

            
        
        
        
        
        
        
        
        
     
        
        
        


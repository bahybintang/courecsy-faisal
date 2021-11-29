import numpy as np
import string
import pandas as pd

import re
import nltk
from nltk.corpus import stopwords
from nltk.tokenize import word_tokenize
from pipeline import lower, remove_punctuation, remove_stopwords, preprocessing, stem_text, lemmatize_text
from sklearn.metrics.pairwise import cosine_similarity

# function untuk membuat bag of words
def createVocab(docs):
    vocab = {} #menyimpan bag of words
    for doc in docs:
        """print("before preprocessed")
        print(doc)"""
        #doc = preprocessing(doc)
        """print("after preprocessed")
        print(doc)
        print()"""
        words = word_tokenize(doc)
        for word in words:
            if(word in vocab.keys()):
                vocab[word] = vocab[word] +1
            else:
                vocab[word] =1 
    return vocab


def TF_IDF_Cosim(docs):
    vocab = createVocab(docs)
    """print("bag of word dan jumlah katanya / TF normal :")
    print(vocab)
    print()"""

    #Compute document term matrix as well idf for each term 
    docsTFMat = np.zeros((len(docs),len(vocab)))

    #docsIdfMat = np.zeros((len(vocab),len(docs)))

    docTermFreq = pd.DataFrame(docsTFMat ,columns=sorted(vocab.keys())) #Ini adalah TF

    docCount=0
    for doc in docs:

        words = word_tokenize(doc)
        for word in words:
            if(word in vocab.keys()):
                docTermFreq[word][docCount] = docTermFreq[word][docCount]+1
            
        docCount = docCount+1

    """print("TF :")    
    print(docTermFreq)
    print()"""

    #Computed idf for each word in vocab
    idfDict={}

    for column in docTermFreq.columns:
        idfDict[column]= np.log((len(docs))/((docTermFreq[column] != 0).sum()))
    """print("nilai IDF tiap term :")    
    print(docTermFreq)
    print()"""

        
    #inisialisasi matriks TF-IDF
    docsTfIdfMat = np.zeros((len(docs),len(vocab)))
    docTfIdf = pd.DataFrame(docsTfIdfMat ,columns=sorted(vocab.keys()))

    #menghitung TF-IDF
    docCount = 0
    for doc in docs:
        for key in idfDict.keys():
            docTfIdf[key][docCount] = docTermFreq[key][docCount] * (idfDict[key]+1)
        docCount = docCount +1 
    """print("Nilai TF-IDF")    
    print(docTfIdf)
    print()"""

    """
    #compute cosine similarity
    csim = cosine_similarity(docTfIdf.values.tolist(),docTfIdf.values.tolist())

    hasil = []
    for x in csim[0]:
        hasil.append(x)

    hasil.pop(0) #hilangkan nilai cosine similarity untuk query sebab tidak dibutuhkan
    """

    hasil = []
    query = docTfIdf.iloc[0,:].to_list()
    for i in range(1, len(docs)):
        video = docTfIdf.iloc[i,:].to_list()
        similarity = np.dot(query, video) / (np.sqrt(np.sum(np.power(query,2))) * np.sqrt(np.sum(np.power(video,2))))
        hasil.append(similarity)

    return hasil
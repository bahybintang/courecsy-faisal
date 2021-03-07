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
        #doc= doc.translate(str.maketrans('', '', string.punctuation))
        #words= word_tokenize(doc.lower())
        doc = preprocessing(doc)
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

def coSim(docs):
    vocab = createVocab(docs)
    """print("bag of word dan jumlah katanya / TF normal :")
    print(vocab)
    print()"""

    #Compute document term matrix as well idf for each term 
    docsTFMat = np.zeros((len(docs),len(vocab)))

    docsIdfMat = np.zeros((len(vocab),len(docs)))

    docTermFreq = pd.DataFrame(docsTFMat ,columns=sorted(vocab.keys())) #Ini adalah TF

    docCount=0
    for doc in docs:
        #doc= doc.translate(str.maketrans('', '', string.punctuation))
        #words= word_tokenize(doc.lower())
        doc = preprocessing(doc)
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
        idfDict[column]= np.log((len(docs) +1 )/(1+ (docTermFreq[column] != 0).sum()))+1
    """print("nilai IDF tiap term :")    
    print(docTermFreq)
    print()"""

        
    #compute tf.idf matrix
    docsTfIdfMat = np.zeros((len(docs),len(vocab)))
    docTfIdfDf = pd.DataFrame(docsTfIdfMat ,columns=sorted(vocab.keys()))

    docCount = 0
    for doc in docs:
        for key in idfDict.keys():
            docTfIdfDf[key][docCount] = docTermFreq[key][docCount] * idfDict[key]
        docCount = docCount +1 
    """print("Nilai TF-IDF")    
    print(docTfIdfDf)
    print()"""

    #compute cosine similarity
    csim = cosine_similarity(docTfIdfDf.values.tolist(),docTfIdfDf.values.tolist())
    """print("Tabel perbandingan cosine similarity :")
    print(csim)
    print()"""

    hasil = []
    for x in csim[0]:
        hasil.append(x)
    """print("Hasil cosine similairty (masih ada query):")
    print(hasil)
    print()"""

    hasil.pop(0) #hilangkan nilai cosine similarity untuk query sebab tidak dibutuhkan
    """print("Hasil cosine similairty tanpa query:")
    print(hasil)
    print()"""

    return hasil
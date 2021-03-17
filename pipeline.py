import re
import nltk
from nltk.corpus import stopwords
from nltk.stem import PorterStemmer
from nltk.stem import WordNetLemmatizer

stp = stopwords.words('english')
stm = PorterStemmer()
lmtz = WordNetLemmatizer()


def lower(text):
    return text.lower()


def remove_punctuation(text):
    # hilangkan url
    text = re.sub(r'(?:(?:https?|ftp):\/\/)?[\w/\-?=%.]+\.[\w/\-&?=%.]+', ' ', text)

    # hilangkan tag html
    text = re.sub(r'<.+?>', ' ', text)

    # hialangkan simbol
    text = re.sub(r'[^\w\s+]', ' ', text)
    return text


def remove_stopwords(text):
    text = ' '.join([word for word in text.split() if word not in stp])
    return text

def preprocessing(text):
    text = lower(text)
    text = remove_punctuation(text)
    text = remove_stopwords(text)
    return text

def stem_text(text):
    text = ' '.join([stm.stem(word) for word in text.split()])
    return text

def preproc_stem(text):
    text = lower(text)
    text = remove_punctuation(text)
    text = remove_stopwords(text)
    text = stem_text(text)
    return text

def lemmatize_text(text):
    text = ' '.join([lmtz.lemmatize(word) for word in text.split()])
    return text
from flask import Flask, render_template, request, redirect, url_for
from wtforms import Form, StringField, TextAreaField, PasswordField, validators
import requests
from isodate import parse_duration, parse_date
import sys
import logging

from pipeline import lower, remove_punctuation, remove_stopwords, preprocessing, stem_text, lemmatize_text
from csim import createVocab, coSim
from dataset import Silabus

app = Flask(__name__)

app.logger.addHandler(logging.StreamHandler(sys.stdout))
app.logger.setLevel(logging.ERROR)

data_silabus = Silabus()

@app.route('/', methods=["GET", "POST"])
def dashboard():
    return render_template('dashboard.html', data_silabus=data_silabus)

"""@app.route('/dashboard')
def dashboard():
    return render_template('dashboard.html', data_silabus=data_silabus)"""


class ArticleForm(Form):
    title = StringField('Title', [validators.Length(min=1, max=200)])
    body = TextAreaField('Body (tulis apabila ada silabus yang ingin ditambahkan)')


@app.route('/edit_silabus/<string:id>', methods=['GET', 'POST'])
def edit_article(id):
    article = data_silabus[int(id)-1]
    # Get form
    form = ArticleForm(request.form)

    # Populate article form fields
    form.title.data = article['title']
    #form.body.data = article['silabus']

    #####################################################################################
    search_url = "https://www.googleapis.com/youtube/v3/search"
    video_url = "https://www.googleapis.com/youtube/v3/videos"

    videos = [] ## menyimpan dictionary dari data video

    if request.method == 'POST' and form.validate():
        title = request.form['title']
        box = request.form.getlist('silabus')
        body = request.form['body']
        durasi = request.form['durasi']
        tanggal = request.form['tanggal']
        
        box_to_text = ' '.join(box)
        query = title + " " + box_to_text + " " + body

        print(type(body))
        print(body)
        print(box)
        print("query yang diinput : " + query)

        if tanggal == "":
            search_params = {
                "key" : "AIzaSyCCmLV6KrNBhsQmpc5CVdsRCvKfrKxayuA",
                "q" : preprocessing(query),
                "part" : "snippet",
                "maxResults" : 50,
                "type" : "video",
                "videoDuration" : durasi
            }
        else :
            search_params = {
                "key" : "AIzaSyCCmLV6KrNBhsQmpc5CVdsRCvKfrKxayuA",
                "q" : preprocessing(query),
                "part" : "snippet",
                "maxResults" : 50,
                "type" : "video",
                "videoDuration" : durasi,
                "publishedBefore" : tanggal + "T00:00:00Z"
            }
        
        print(search_params["q"])

        r = requests.get(search_url, params=search_params)

        #print(r.text) ##tampilkan semua atribut video
        results = r.json()["items"]

        video_ids = []

        for result in results:
            video_ids.append(result["id"]["videoId"])

        video_params = {
            "key" : "AIzaSyCCmLV6KrNBhsQmpc5CVdsRCvKfrKxayuA",
            "id" : ",".join(video_ids),
            "part" : "snippet, contentDetails",
            "maxResults" : 50
        }
        
        r = requests.get(video_url, params=video_params)
        results = r.json()["items"]

        for result in results:
            video_data = {
                "id" : result["id"],
                "url" : f'https://www.youtube.com/watch?v={result["id"]}',
                "thumbnail" : result["snippet"]["thumbnails"]["high"]["url"],
                "duration" : int(parse_duration(result["contentDetails"]["duration"]).total_seconds() // 60),
                "title" : result["snippet"]["title"],
                "description" : result["snippet"]["title"] + " " + result["snippet"]["description"],
                "published" : parse_date(result["snippet"]["publishedAt"])
            }
            # Urutkan berdasarkan nilai cosine similaritynya
            videos.append(video_data)
        
        #print(search_params["q"])
        #rint(videos[0]["description"])
        #print(videos)

        docs = [] #list of query dan deskripsi semua video
        docs.append(search_params["q"])

        for video in videos:
            docs.append(preprocessing(video["description"]))

        print(docs)
        
        cosine_similairty = coSim(docs)

        i=0
        for video in videos:
            video["similarity"] = cosine_similairty[i]
            i += 1

        #return render_template('index.html', videos=sorted(videos, key = lambda i: i["similarity"], reverse=True))
        return render_template('rekomendasi.html', videos=sorted(videos, key = lambda i: i["similarity"], reverse=True), title=title)
    ######################################################################################
    return render_template('edit_silabus.html', form=form, article=article)

@app.route('/rekomendasi')
def rekomendasi():
    return render_template('rekomendasi.html')

if __name__ == "__main__":
    app.run(debug = True)

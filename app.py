from flask import Flask, render_template, request, redirect, flash
from wtforms import Form, StringField, TextAreaField, validators
import requests
from isodate import parse_duration, parse_date #converting the ISO_8601 format

from pipeline import preprocessing
from csim import TF_IDF_Cosim
from dataset import Silabus

app = Flask(__name__)

data_silabus = Silabus()


@app.route("/", methods=["GET", "POST"])
def dashboard():
    return render_template("dashboard.html", data_silabus=data_silabus)

class ArticleForm(Form):
    title = StringField("Title", [validators.Length(min=1, max=200)])
    body = TextAreaField("Body (tulis apabila ada silabus yang ingin ditambahkan)")


@app.route("/edit_silabus/<string:id>", methods=["GET", "POST"])
def edit_article(id):
    article = data_silabus[int(id) - 1]
    # Get form
    form = ArticleForm(request.form)

    # Populate article form fields
    form.title.data = article["title"]
    # form.body.data = article['silabus']

    #####################################################################################
    videos_data = []  ## untuk menyimpan data video yang didapat dari YouTube

    if request.method == "POST" and form.validate():
        title = request.form["title"]
        box = request.form.getlist("silabus")
        body = request.form["body"]
        durasi = request.form["durasi"]
        tanggal = request.form["tanggal"]

        if len(box) == 0 and body == "":
            flash(
                "Form checkbox silabus dan body masih kosong, mohon isi salah satu atau keduanya",
                "danger",
            )
            return redirect(id)

        box_to_text = " ".join(box)
        query = title + " " + box_to_text + " " + body
        """query = "numerical method interpolation taylor series definitions divided differences newton lagrange errors polynomial chebyshev hermite spline """

        print("Input mata kuliah dan silabus : \n" + query + "\n")
        print(
            "Hasil pemrosesan teks pada input mata kuliah dan silabus : \n"
            + preprocessing(query)
        )

        if tanggal == "":
            search_params = {
                "key": "AIzaSyCCmLV6KrNBhsQmpc5CVdsRCvKfrKxayuA",
                "q": preprocessing(query),
                "part": "snippet",
                "maxResults": 50,
                "type": "video",
                "videoDuration": durasi,
            }
        else:
            search_params = {
                "key": "AIzaSyCCmLV6KrNBhsQmpc5CVdsRCvKfrKxayuA",
                "q": preprocessing(query),
                "part": "snippet",
                "maxResults": 50,
                "type": "video",
                "videoDuration": durasi,
                "publishedBefore": tanggal + "T00:00:00Z",
            }

        search_url = "https://www.googleapis.com/youtube/v3/search"
        r = requests.get(search_url, params=search_params)

        results = r.json()["items"]

        video_ids = []

        for result in results:
            video_ids.append(result["id"]["videoId"])

        video_params = {
            "key": "AIzaSyCCmLV6KrNBhsQmpc5CVdsRCvKfrKxayuA",
            "id": ",".join(video_ids),
            "part": "snippet, contentDetails",
            "maxResults": 10,
        }

        video_url = "https://www.googleapis.com/youtube/v3/videos"
        r = requests.get(video_url, params=video_params)
        results = r.json()["items"]

        judul_video = []

        for result in results:
            video_data = {
                "id": result["id"],
                "url": f'https://www.youtube.com/watch?v={result["id"]}',
                "thumbnail": result["snippet"]["thumbnails"]["high"]["url"],
                "duration": int(
                    parse_duration(result["contentDetails"]["duration"]).total_seconds()
                    // 60
                ),
                "title": result["snippet"]["title"],
                "description": result["snippet"]["title"]
                + " "
                + result["snippet"]["description"],
                "published": parse_date(result["snippet"]["publishedAt"]),
            }
            # Urutkan berdasarkan nilai cosine similaritynya
            videos_data.append(video_data)
            judul_video.append(video_data["title"])

        print(judul_video)

        docs = []  # list dari teks query dan deskripsi semua video
        docs.append(search_params["q"])

        for video in videos_data:
            docs.append(preprocessing(video["description"]))

        cosine_similairty = TF_IDF_Cosim(docs)

        i = 0
        for video in videos_data:
            video["similarity"] = "{:.6f}".format(cosine_similairty[i])
            i += 1

        return render_template(
            "rekomendasi.html",
            videos_data=sorted(
                videos_data, key=lambda i: i["similarity"], reverse=True
            ),
            title=title,
        )

    return render_template("edit_silabus.html", form=form, article=article)


# @app.route("/rekomendasi")
# def rekomendasi():
#     return render_template("rekomendasi.html")


if __name__ == "__main__":
    app.secret_key = "secret123"
    app.run(debug=True)

FROM python:3.8-slim-buster
ENV FLASK_APP app
ENV FLASK_ENV production
WORKDIR /app
COPY . .
RUN pip3 install -r requirements.txt
RUN python3 -m nltk.downloader stopwords
CMD [ "python3", "-m" , "flask", "run", "--host=0.0.0.0"]

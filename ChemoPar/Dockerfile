FROM python:3.11.8-slim-bookworm

ENV PYTHONBUFFERED=1
ENV SECRET_KEY='78g981h3fi2n982fehoi'

ENV PORT 8080

WORKDIR /app

COPY . /app/

RUN pip install --upgrade pip
RUN pip install -r requirements.txt

RUN python manage.py collectstatic --noinput

CMD gunicorn ChemoPar.wsgi:application --bind 0.0.0.0:"${PORT}"

EXPOSE ${PORT} 

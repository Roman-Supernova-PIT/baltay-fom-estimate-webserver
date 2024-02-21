# build with
#   docker build -t baltay-fom-webserver .
#
# run with
#   docker run -d --name baltay-fom -p 8080:8080 baltay-fom-webserver
#
# To bind-mount the source directory for testing purposes, after -d:
#    --mount type=bind,source=.,target=/code

FROM docker.io/alpine:3.19.1

RUN apk update \
    && apk add python3 py3-pip \
    && apk cache clean \
    && rm -rf /var/cache/apk/*

RUN rm /usr/lib/python3.11/EXTERNALLY-MANAGED

RUN pip install gunicorn flask \
    && rm -rf /.cache/pip

RUN mkdir /code
RUN mkdir /code/static
RUN mkdir /code/templates
WORKDIR /code
COPY webservice.py webservice.py
COPY static/*.js /code/static/
COPY rkwebutil/rkwebutil.js /code/static/rkwebutil.js
COPY templates/*.html /code/templates/

CMD [ "gunicorn", "-w", "4", "-b", "0.0.0.0:8080", "--timeout", "0", \
      "webservice:create_app()" ]

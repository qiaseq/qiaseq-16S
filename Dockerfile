FROM python:2.7.15-alpine
MAINTAINER Raghavendra Padmanabhan <raghavendra.padmanabhan@qiagen.com>

RUN apk update && apk add alpine-sdk pigz

RUN pip install luigi edlib

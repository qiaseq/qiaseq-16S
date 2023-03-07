FROM python:2.7.15-alpine

LABEL author="Raghavendra Padmanabhan"
LABEL maintainer="QIAGEN Biox" email="LSBioXRequests@qiagen.com"

RUN apk update && apk add alpine-sdk pigz

RUN pip install luigi edlib

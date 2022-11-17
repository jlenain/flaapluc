#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Time-stamp: "2022-11-17 12:13:36 jlenain"

"""
Alert producer for FLaapLUC using the kafka protocol and AVRO alert format
"""

import yaml
import os

from confluent_kafka import Producer
from flaapluc import avroUtils

import flaapluc


class AlertProducer:

    def __init__(self, conf_path=None):
        try:
            with open(conf_path) as f:
                self.conf = yaml.load(f, Loader=yaml.FullLoader)
        except FileNotFoundError as e:
            print(e)
            raise

        self.topic = self.conf['topic']
        self.prod_conf = {'bootstrap.servers': self.conf['server'],
                          'group.id': self.conf['group']}
        self.schema = f"{os.path.dirname(flaapluc.__file__)}/../schemas/{self.conf['schema']}"

    def sendAlert(self, alert=None):
        p = Producer(self.prod_conf)
        avro_alert = avroUtils.encode_into_avro(alert, self.schema)
        p.produce(self.topic, avro_alert)

        p.flush()

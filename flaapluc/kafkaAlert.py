#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Time-stamp: "2022-11-16 12:06:23 jlenain"

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

    def sendAlert(self, alerts=None):
        p = Producer(self.prod_conf)

        for alert in alerts:
            print(alert)
            avro_alert = avroUtils.encode_into_avro(alert, self.schema)
            p.produce(self.topic, avro_alert)

        p.flush()


# Test AlertProducer
conf_path = f"{os.path.dirname(flaapluc.__file__)}/../conf/flaapluc-kafka-broker-configuration.yml"

a = AlertProducer(conf_path=conf_path)
alerts = []
for i in range(50):
    alert = dict()
    alert['alert'] = dict()
    msg = f'Hi from CC-IN2P3 ! Message number {i}.'
    alert['alert']['comment'] = msg
    alerts.append(alert)

a.sendAlert(alerts)

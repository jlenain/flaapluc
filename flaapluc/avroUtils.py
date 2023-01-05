#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Time-stamp: "2022-11-17 12:16:18 jlenain"

import io
import json

import fastavro


# From https://github.com/astrolabsoftware/fink-client/blob/master/fink_client/avroUtils.py
def encode_into_avro(alert: dict, schema_file: str) -> str:
    """Encode a dict record into avro bytes
    Parameters
    ----------
    alert: dict
        A Dictionary of alert data
    schema_file: str
        Path of avro schema file
    Returns
    ----------
    value: str
        a bytes string with avro encoded alert data
    """
    with open(schema_file) as f:
        schema = json.load(f)

    parsed_schema = fastavro.parse_schema(schema)
    b = io.BytesIO()
    fastavro.schemaless_writer(b, parsed_schema, alert)

    return b.getvalue()

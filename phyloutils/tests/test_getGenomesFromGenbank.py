#!/usr/bin/env python

import os
import argparse
import logging
from unittest import TestCase
from phyloutils.subcommands import getGenomesFromGenbank


def download_genbank(out_dir="test_output"):
    input_file = "./data/assembly_result.txt"
    out_dir = out_dir
    parser = argparse.ArgumentParser()
    parser = getGenomesFromGenbank.add_arguments(parser)
    args = parser.parse_args(['-t', 'ID', '-o', out_dir, input_file])
    getGenomesFromGenbank.command(args)


class TestGetGenomesFromGenbank(TestCase):
    def test_download_genbank(self):
        out_dir = "test_output"
        download_genbank(out_dir)
        abs_path = os.path.abspath(out_dir)
        for folder in os.listdir(out_dir):
            if os.path.isdir(folder):
                self.assertIsNotNone(os.listdir(os.path.join(abs_path, folder)))
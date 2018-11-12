import os
import xml.etree.cElementTree as et
import re
import pandas as pd


def import_trials():
    """
    Given a set of clinical trial files, read them in.
    """

    # read all trial source files (*.xml) in trial data sub-directories
    # total num of files must be 241006
    trial_files = []
    for root, dirs, files in os.walk('../../data/clinicaltrials_txt/'):
        if not re.match(r"\d+", root.split(os.sep)[-1]):
            continue
        for file in files:
            if not file.endswith('txt'):
                continue
            trial_files.append(os.path.join(root, file))
    print('INFO', '{} trials found'.format(len(trial_files)))

    trials = []
    for i, file in enumerate(trial_files):
        # print('INFO', 'parsing a trial file {}'.format(file))
        with open(file, 'r') as myfile:
            data = myfile.read().replace('\n', '')
        trials.append(data.lower())

    return trials


def getvalueofnode(node):
    """ return node text or None """
    return node.text if node is not None else None


def import_topics():
    """
    Read in a set of patient profiles (topics)
    """

    topics = []
    parsedXML = et.parse("../../data/topics2017.xml")
    for node in parsedXML.getroot():
        disease = getvalueofnode(node.find('disease')).lower()
        gene = getvalueofnode(node.find('gene')).lower()
        demographic = getvalueofnode(node.find('demographic')).lower()
        other = getvalueofnode(node.find('other')).lower()
        topics.append((disease, gene, demographic, other))

    return topics


def count_usages(trials, topics):
    """
    Given a set of trials and topics, return a dataset giving counts of times
    diseases, genes, demographics, or other value occurred in each pairing of
    trials and topics.
    """
    trials_rows = []
    for trial in trials:
        counts = []
        for topic in topics:
            disease = topic[0]
            count = 0
            if disease in trial:
                count = count + 1
            counts.append(count)
        trials_rows.append(counts)
    df = pd.DataFrame(trials_rows)
    return(df)

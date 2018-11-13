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
            data = myfile.read()
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
    Given a set of trials and topics, return a dataset giving score based
    number of times of times
    diseases, genes, demographics, or other value occurred in each pairing of
    trials and topics.
    """
    trials_rows = []
    for trial in trials:
        scores = []
        for topic in topics:
            disease = topic[0]
            genes = topic[1].split(', ')
            demographic = topics[3]
            other = topics[4]
            score = trial.count(disease)
            match = True
            if score == 0:
                match = False
            score_before = score
            for gene in genes:
                score = score + trial.count(gene)
            if score == score_before:
                match = False
            if match is False:
                score = 0
            scores.append(score)
        trials_rows.append(scores)
    df = pd.DataFrame(trials_rows)
    return(df)


def print_trial(trials, row_index):
    """
    Save trial to text file.
    """
    file_name = "trial" + str(row_index) + '.txt'
    with open(file_name, "w") as text_file:
        text_file.write(trials[row_index])

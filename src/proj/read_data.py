import os
import xml.etree.cElementTree as et
import re
import pandas as pd
from Bio import Entrez
Entrez.email = "lgw2@uw.edu"


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
    trial_names = []
    for i, file in enumerate(trial_files):
        with open(file, 'r') as myfile:
            data = myfile.read()
        trials.append(data.lower())
        trial_names.append(file[40:-4])

    df = pd.DataFrame({
                       'trial_name': trial_names,
                       'trial': trials
                      })

    return df


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


def get_trial_gender(trial):
    trial_gender = trial.split(
        'eligibility:\ngender: '
    )[1].split('\nage:')[0]
    if trial_gender == 'all':
        trial_gender = ['female', 'male']
    else:
        trial_gender = [trial_gender]
    return(trial_gender)


def get_trial_age_to(trial):
    trial_age = trial.split(
        '\nage: '
    )[1].split('\n')[0]
    if trial_age == '20-40years':
        trial_age = '20 years to 40 years'
    age_to = trial_age.split(' to ')[1].split(' year')[0]
    if (age_to == 'n/a') or (age_to == 'any'):
        age_to = 150
    elif ('month' in age_to) or ('day' in age_to) or\
            ('week' in age_to) or ('hour' in age_to) or\
            ('minute' in age_to):
        age_to = 2
    else:
        age_to = int(age_to)
    return(age_to)


def get_trial_age_from(trial):
    trial_age = trial.split(
        '\nage: '
    )[1].split('\n')[0]
    age_from = trial_age.split(' to ')[0].split(' year')[0].split('-')[0]
    if (age_from == 'n/a') or ('week' in age_from) or\
            ('month' in age_from) or ('day' in age_from) or\
            (age_from == 'any'):
        age_from = 0
    else:
        age_from = int(age_from)
    return(age_from)


def get_trial_exclusion(trial):
    trial_sections = trial.split('\n\n        exclusion criteria:\n\n')
    if len(trial_sections) == 2:
        trial_exclusion = trial_sections[1]
    else:
        trial_exclusion = ""
    return(trial_exclusion)


def match_gender(gender, trial_gender, match):
    if gender not in trial_gender:
        return(False)
    else:
        return(match)


def match_age(age, age_from, age_to, match):
    if (age < age_from) or (age > age_to):
        return(False)
    else:
        return(match)


def get_disease(topic):
    return(topic[0])


def get_genes(topic):
    return(topic[1].split(', '))


def get_demographic(topic):
    return(int(topic[2].split('-')[0]), topic[2].split(' ')[1])


def get_other_conditions(topic):
    return(topic[3].split(', '))


def match_genes(genes, trial, score, match):
    score_before = score
    for gene in genes:
        score = score + trial.count(gene)
    if score == score_before:
        match = False
    return(score, match)


def match_other_conds(other_conditions, trial_exclusion, match):
    if other_conditions != ['none']:
        for cond in other_conditions:
            if cond in trial_exclusion:
                match = False
    return(match)


def compute_baseline_scores(trials, topics):
    """
    Given a set of trials and topics, return a dataset giving score based
    number of times of times
    diseases, genes, demographics, or other value occurred in each pairing of
    trials and topics.
    """
    trials_rows = []
    for trial in trials['trial']:
        scores = []
        trial_gender = get_trial_gender(trial)
        age_from = get_trial_age_from(trial)
        age_to = get_trial_age_to(trial)
        trial_exclusion = get_trial_exclusion(trial)
        for topic in topics:
            disease = get_disease(topic)
            genes = get_genes(topic)
            age, gender = get_demographic(topic)
            other_conditions = get_other_conditions(topic)
            score = trial.count(disease)
            match = True
            match = match_gender(gender, trial_gender, match)
            match = match_age(age, age_from, age_to, match)
            if score == 0:
                match = False
            score, match = match_genes(genes, trial, score, match)
            match = match_other_conds(other_conditions, trial_exclusion,
                                      match)
            if match is False:
                score = 0
            scores.append(score)
        trials_rows.append(scores)
    df = pd.DataFrame(trials_rows)
    df['trial'] = trials['trial_name']
    return(df)


def print_trial(trials, row_index):
    """
    Save trial to text file.
    """
    file_name = "trial" + str(row_index) + '.txt'
    with open(file_name, "w") as text_file:
        text_file.write(trials[row_index])


def import_ground_truth():
    file = "../../data/ground_truth.txt"
    data = pd.read_csv(file, sep=' ', header=None,
                       names=['topic', '1', 'trial', '3'])
    return(data[['topic', 'trial']])


def generate_baseline_data(baseline, ground_truth):
    topics = []
    counts = []
    for topic in range(1, 31):
        count = 0
        for score in (set(baseline[topic-1]) - set([0])):
            trials = [x for x in baseline[baseline[topic-1] == score]['trial']]
            count = count + len(trials)
        topics.append(topic)
        counts.append(count)
    return(pd.DataFrame({
                         'topic': topics,
                         'matches': counts
                         }))


def compute_baseline_evaluation(baseline_data, ground_truth):
    prec5 = baseline_data['matches'].sum()/(5*(baseline_data.shape[0]-1))
    prec10 = baseline_data['matches'].sum()/(10*(baseline_data.shape[0]-1))
    prec15 = baseline_data['matches'].sum()/(15*(baseline_data.shape[0]-1))
    overall_precision = 1
    overall_recall = baseline_data['matches'].sum()/ground_truth.shape[0]
    # recall
    evaluation_scores = {
                         'prec5': prec5,
                         'prec10': prec10,
                         'prec15': prec15,
                         'overall_precision': overall_precision,
                         'overall_recall': overall_recall
                        }
    return(evaluation_scores)


def get_gene_aliases(gene):
    print('Gene is: {}'.format(gene))
    handle = Entrez.esearch(db="gene",
                            term="(" + gene +
                            "[Gene Name]) AND homo sapiens[Organism]")
    data = Entrez.read(handle)
    ids = data['IdList']
    print('ids are: {}'.format(ids))
    aliases = []
    for i in ids:
        handle = Entrez.epost(db="gene", id=i)
        result = Entrez.read(handle)
        webEnv = result["WebEnv"]
        queryKey = result["QueryKey"]
        data = Entrez.esummary(db="gene", webenv=webEnv, query_key=queryKey)
        annotations = Entrez.read(data)
        other_aliases = annotations['Document' +
                                    'SummarySet']['Document' +
                                                  'Summary'][0]['OtherAliases']
        print('other aliases are: {}'.format(other_aliases))
        aliases.append(other_aliases)
    return(aliases)


def expand_genes(topics):
    expanded = []
    for topic in topics:
        for gene in [y.split(' ')[0] for y
                     in [x for x in topic[1].split(', ')]]:
            expanded_gene_names = [gene] + get_gene_aliases(gene)
        expanded.append(expanded_gene_names)
    return(expanded)

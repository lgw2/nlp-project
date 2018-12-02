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


def match_expanded_genes(genes, expanded_gene_names, trial, score, match):
    score_before = score
    for gene in genes:
        for expanded_name in expanded_gene_names[gene]:
            score = score + trial.count(expanded_name)
    if score == score_before:
        match = False
    return(score, match)


def match_expanded_diseases(expanded_diseases, trial):
    score = 0
    for disease in expanded_diseases:
        for expanded_name in expanded_diseases[disease]:
            score = score + trial.count(expanded_name)
    return(score)


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


def compute_expanded_scores(trials, topics, expanded_gene_names,
                            expanded_disease_names):
    """
    Given a set of trials and topics, return a dataset giving score based
    number of times of times EXPANDED
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
            score, match = match_expanded_genes(genes, expanded_gene_names,
                                                trial, score, match)
            match = match_other_conds(other_conditions, trial_exclusion,
                                      match)
            if match is False:
                score = 0
            scores.append(score)
        trials_rows.append(scores)
    df = pd.DataFrame(trials_rows)
    df['trial'] = trials['trial_name']
    return(df)



def print_trial(trials, trial_name):
    """
    Save trial to text file.
    """
    file_name = "trial" + trial_name + '.txt'
    print(file_name)
    with open(file_name, "w") as text_file:
        text_file.write(
              trials[trials['trial_name'] == trial_name]['trial'].values[0]
                       )


def import_ground_truth():
    file = '../../data/clinical_trials.judgments.2017.csv'
    data = pd.read_csv(file)
    print("All data: {}".format(data.shape[0]))
    data = data[data['pm_rel_desc'] != 'Not PM']
#    print('After removing not PM: {}'.format(data.shape[0]))
    data = data[data['demographics_desc'] == 'Matches']
    print('after removing demographic mismatches: {}'.format(data.shape[0]))
    data = data[(data['other_desc'] == 'Not Discussed') |
                         (data['other_desc'] == 'Matches')]
    print('after removing other mismatches: {}'.format(data.shape[0]))
    data = data[data['disease_desc'] != 'Not Disease']
    print('after removing disease mismatches: {}'.format(data.shape[0]))
    gene_match = []
    for row_idx in range(data.shape[0]):
        row = data.iloc[row_idx]
        match = False
        for gene in ['1', '2', '3']:
            col_name = 'gene' + gene + '_annotation_desc'
            if row[col_name] in ['Exact',
                                 'Missing Variant',
                                 'Different Variant']:
                match = True
        gene_match.append(match)
    data['gene_match'] = gene_match
    data = data[data['gene_match']]
    print('after removing gene mismatches: {}'.format(data.shape[0]))
    return(data)



def generate_ordered_trials(baseline, ground_truth, topic):
    trials_list = []
    scores = set(baseline[topic-1])-set([0])
    scores = [x for x in scores]
    scores.sort(reverse=True)
    for score in scores:
        print(score)
        trials_at_score = [x for x in baseline[baseline[topic-1] ==\
                  score]['trial']]
        print(trials_at_score)
        trials_list.extend(trials_at_score)
        if len(trials_list) >= 15:
            break
    trials_list = trials_list[:15]
    num_missing_trials = 15 - len(trials_list)
    to_add = baseline.sample(
                    frac=num_missing_trials/baseline.shape[0],
                    random_state=1
                   )['trial'].values
    trials_list.extend(to_add)
    trials_list = trials_list[:15]
    return(trials_list)


def find_num_correct_trials(trials_list, filtered_ground_truth):
    return(len(set(trials_list).intersection(set(filtered_ground_truth))))


def generate_result_data(baseline, ground_truth):
    topics = []
    true_positives_5 = []
    false_positives_5 = []
    true_positives_10 = []
    false_positives_10 = []
    true_positives_15 = []
    false_positives_15 = []
    for topic in range(1, 31):
        topics.append(topic)
        print(topic)
        filtered_ground_truth =\
               ground_truth[ground_truth['trec_topic_number'] ==
               topic]['trec_doc_id']
        trials = generate_ordered_trials(baseline, ground_truth, topic)
        first_5 = trials[:5]
        first_10 = trials[:10]
        first_15 = trials[:15]
        tp5 = find_num_correct_trials(first_5, filtered_ground_truth)
        fp5 = 5 - tp5
        tp10 = find_num_correct_trials(first_10, filtered_ground_truth)
        fp10 = 10 - tp10
        tp15 = find_num_correct_trials(first_15, filtered_ground_truth)
        fp15 = 15 - tp15
        true_positives_5.append(tp5)
        false_positives_5.append(fp5)
        true_positives_10.append(tp10)
        false_positives_10.append(fp10)
        true_positives_15.append(tp15)
        false_positives_15.append(fp15)
    return(pd.DataFrame({
                         'topic': topics,
                         'true_positives_5': true_positives_5,
                         'false_positives_5': false_positives_5,
                         'true_positives_10': true_positives_10,
                         'false_positives_10': false_positives_10,
                         'true_positives_15': true_positives_15,
                         'false_positives_15': false_positives_15
                         }))


def compute_result_evaluation(baseline_data, ground_truth):
    num_topics = 29
    prec5 = baseline_data['true_positives_5'].sum()/(5*num_topics)
    prec10 = baseline_data['true_positives_10'].sum()/(10*num_topics)
    prec15 = baseline_data['true_positives_15'].sum()/(15*num_topics)
    evaluation_scores = {
                         'prec5': prec5,
                         'prec10': prec10,
                         'prec15': prec15,
                        }
    return(evaluation_scores)


def get_gene_aliases(gene):
    handle = Entrez.esearch(db="gene",
                            term="(" + gene +
                            "[Gene Name]) AND homo sapiens[Organism]")
    data = Entrez.read(handle)
    ids = data['IdList']
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
        other_aliases = other_aliases.split(', ')
        aliases.extend(other_aliases)
    return(aliases)


def expand_genes(topics):
    expanded = {}
    for topic in topics:
        for gene in topic[1].split(', '):
            simple_gene_name = gene.split(' ')[0]
            expanded_gene_names = [simple_gene_name] + get_gene_aliases(gene)
            expanded[gene] = expanded_gene_names
    return(expanded)


def import_doid():
    doid = pd.read_csv('../../data/DOID.csv')
    return(doid)

def expand_diseases(topics, doid):
    expanded = {}
    for topic in topics:
        disease = topic[0]
        print(disease)
        synonyms = doid[doid['Preferred Label'] == disease]['Synonyms']
        if synonyms.shape[0] > 0:
            if not synonyms.isnull().values[0]:
                syn_list = synonyms.values[0]
                syn_list = syn_list.split('|')
                expanded[disease] = [disease] + syn_list
            else:
                expanded[disease] = [disease]
        else:
            expanded[disease] = [disease]
    return(expanded)
